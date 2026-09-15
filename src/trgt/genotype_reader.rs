use crate::utils::Result;
use rust_htslib::{
    bam::{self, Read, Record},
    tpool::ThreadPool,
};
use std::sync::Arc;

pub trait BamRead {
    fn fetch(&mut self, region: (&str, u32, u32)) -> Result<()>;
    fn read(&mut self, record: &mut Record) -> Option<Result<()>>;
}

pub enum BamReader {
    Single(SingleBamReader),
    Multi(MultiBamReader),
}

impl BamReader {
    pub fn single(reader: bam::IndexedReader) -> Self {
        Self::Single(SingleBamReader::new(reader))
    }

    pub fn multi(sources: Vec<BamSourceReader>) -> Self {
        Self::Multi(MultiBamReader::new(sources))
    }

    pub fn multi_with_tpool(sources: Vec<BamSourceReader>, thread_pool: ThreadPool) -> Self {
        Self::Multi(MultiBamReader::with_tpool(sources, thread_pool))
    }

    pub fn fetch(&mut self, region: (&str, u32, u32)) -> Result<()> {
        BamRead::fetch(self, region)
    }

    pub fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        BamRead::read(self, record)
    }
}

impl BamRead for BamReader {
    fn fetch(&mut self, region: (&str, u32, u32)) -> Result<()> {
        match self {
            Self::Single(reader) => reader.fetch(region),
            Self::Multi(reader) => reader.fetch(region),
        }
    }

    fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        match self {
            Self::Single(reader) => reader.read(record),
            Self::Multi(reader) => reader.read(record),
        }
    }
}

pub struct BamSourceReader {
    source_label: Arc<str>,
    reader: bam::IndexedReader,
}

impl BamSourceReader {
    pub fn new(source_label: impl Into<Arc<str>>, reader: bam::IndexedReader) -> Self {
        Self {
            source_label: source_label.into(),
            reader,
        }
    }
}

pub struct SingleBamReader {
    reader: bam::IndexedReader,
}

impl SingleBamReader {
    pub fn new(reader: bam::IndexedReader) -> Self {
        Self { reader }
    }
}

impl BamRead for SingleBamReader {
    fn fetch(&mut self, region: (&str, u32, u32)) -> Result<()> {
        self.reader.fetch(region).map_err(|e| e.to_string())
    }

    fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        self.reader
            .read(record)
            .map(|result| result.map_err(|e| e.to_string()))
    }
}

pub struct MultiBamReader {
    sources: Vec<BufferedSourceReader>,
    // Source readers are declared before the pool so they are dropped first.
    _tpool: Option<ThreadPool>,
}

impl MultiBamReader {
    pub fn new(sources: Vec<BamSourceReader>) -> Self {
        Self::new_inner(sources, None)
    }

    pub fn with_tpool(sources: Vec<BamSourceReader>, tpool: ThreadPool) -> Self {
        Self::new_inner(sources, Some(tpool))
    }

    fn new_inner(sources: Vec<BamSourceReader>, tpool: Option<ThreadPool>) -> Self {
        let sources = sources.into_iter().map(BufferedSourceReader::new).collect();
        Self {
            sources,
            _tpool: tpool,
        }
    }
}

impl BamRead for MultiBamReader {
    fn fetch(&mut self, region: (&str, u32, u32)) -> Result<()> {
        let (contig, start, end) = region;
        for source in &mut self.sources {
            if let Err(e) = source.reader.fetch(region) {
                return Err(format!(
                    "BAM fetch error for source {} querying {contig}:{start}-{end}: {e}",
                    source.source_label.as_ref(),
                ));
            }
            source.reset();
        }
        Ok(())
    }

    fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        let source_index = match self.next_source_index() {
            Ok(Some(source_index)) => source_index,
            Ok(None) => return None,
            Err(e) => return Some(Err(e)),
        };
        let source = &mut self.sources[source_index];
        source
            .current_record_sort_key
            .take()
            .expect("next_source_index returns only sources with a current record");
        std::mem::swap(record, &mut source.record);
        Some(Ok(()))
    }
}

impl MultiBamReader {
    fn next_source_index(&mut self) -> Result<Option<usize>> {
        let mut next = None;

        // For the current reads + fail reads option, a linear scan is cheaper than
        // maintaining a heap. If we ever support many BAM inputs, revisit this with a
        // BinaryHeap keyed by current_sort_key().
        for (index, source) in self.sources.iter_mut().enumerate() {
            if source.current_record_sort_key.is_none() && !source.is_exhausted {
                source.fill_current_record()?;
            }

            let Some(key) = source.current_sort_key(index) else {
                continue;
            };

            match next {
                Some((_, next_key)) if next_key <= key => {}
                _ => next = Some((index, key)),
            }
        }

        Ok(next.map(|(index, _)| index))
    }
}

struct BufferedSourceReader {
    source_label: Arc<str>,
    reader: bam::IndexedReader,
    record: Record,
    current_record_sort_key: Option<(i32, i64, u64)>,
    is_exhausted: bool,
    next_record_order: u64,
}

impl BufferedSourceReader {
    fn new(source: BamSourceReader) -> Self {
        Self {
            source_label: source.source_label,
            reader: source.reader,
            record: Record::new(),
            current_record_sort_key: None,
            is_exhausted: false,
            next_record_order: 0,
        }
    }

    fn reset(&mut self) {
        self.current_record_sort_key = None;
        self.is_exhausted = false;
        self.next_record_order = 0;
    }

    fn fill_current_record(&mut self) -> Result<()> {
        match self.reader.read(&mut self.record) {
            Some(Ok(())) => {
                let order = self.next_record_order;
                self.next_record_order += 1;
                self.current_record_sort_key = Some((self.record.tid(), self.record.pos(), order));
                Ok(())
            }
            Some(Err(e)) => Err(format!(
                "BAM read error for source {}: {}",
                self.source_label.as_ref(),
                e
            )),
            None => {
                self.is_exhausted = true;
                Ok(())
            }
        }
    }

    #[inline]
    fn current_sort_key(&self, source_index: usize) -> Option<(i32, i64, usize, u64)> {
        self.current_record_sort_key
            .map(|(tid, pos, order)| (tid, pos, source_index, order))
    }
}

#[cfg(test)]
mod tests {
    use super::{BamRead, BamReader, BamSourceReader, BufferedSourceReader, MultiBamReader};
    use crate::utils::test_util::create_indexed_bam;
    use rust_htslib::bam::{self, Read, Record};
    use tempfile::NamedTempFile;

    fn source_reader(label: &str, bam_file: &NamedTempFile) -> BamSourceReader {
        let reader = bam::IndexedReader::from_path(bam_file.path()).unwrap();
        BamSourceReader::new(label, reader)
    }

    fn multi_reader(sources: &[(&str, &NamedTempFile)]) -> BamReader {
        BamReader::multi(
            sources
                .iter()
                .map(|(label, bam_file)| source_reader(label, bam_file))
                .collect(),
        )
    }

    fn collect_records(reader: &mut BamReader) -> Vec<(String, i64)> {
        let mut records = Vec::new();
        let mut record = Record::new();
        while let Some(result) = reader.read(&mut record) {
            result.expect("failed to read BAM record");
            records.push((
                std::str::from_utf8(record.qname()).unwrap().to_owned(),
                record.pos(),
            ));
        }
        records
    }

    fn collect_qnames(reader: &mut BamReader) -> Vec<String> {
        let mut qnames = Vec::new();
        let mut record = Record::new();
        while let Some(result) = reader.read(&mut record) {
            result.expect("failed to read BAM record");
            qnames.push(std::str::from_utf8(record.qname()).unwrap().to_owned());
        }
        qnames
    }

    #[test]
    fn single_source_reader_delegates_fetch_and_read() {
        let bam_file = create_indexed_bam(&[
            ("before-region", "chr1", 10),
            ("in-region-1", "chr1", 100),
            ("in-region-2", "chr1", 130),
            ("after-region", "chr1", 300),
        ]);
        let mut direct = bam::IndexedReader::from_path(bam_file.path()).unwrap();
        direct.fetch(("chr1", 50, 200)).unwrap();
        let mut expected = Vec::new();
        let mut record = Record::new();
        while let Some(result) = direct.read(&mut record) {
            result.unwrap();
            expected.push((
                std::str::from_utf8(record.qname()).unwrap().to_owned(),
                record.pos(),
            ));
        }

        let indexed_reader = bam::IndexedReader::from_path(bam_file.path()).unwrap();
        let mut bam_reader = BamReader::single(indexed_reader);
        bam_reader.fetch(("chr1", 50, 200)).unwrap();

        assert_eq!(collect_records(&mut bam_reader), expected);
        assert_eq!(
            expected,
            vec![
                ("in-region-1".to_string(), 100),
                ("in-region-2".to_string(), 130)
            ]
        );
    }

    #[test]
    fn multi_source_fetch_reads_the_interval_from_every_source() {
        let primary = create_indexed_bam(&[("primary-in-region", "chr1", 100)]);
        let failed = create_indexed_bam(&[("failed-in-region", "chr1", 110)]);
        let mut reader = multi_reader(&[("primary.bam", &primary), ("failed.bam", &failed)]);

        reader.fetch(("chr1", 50, 200)).unwrap();

        assert_eq!(
            collect_qnames(&mut reader),
            vec!["primary-in-region", "failed-in-region"]
        );
    }

    #[test]
    fn multi_source_records_are_emitted_in_virtual_coordinate_order() {
        let source_0 =
            create_indexed_bam(&[("source-0-100", "chr1", 100), ("source-0-120", "chr1", 120)]);
        let source_1 =
            create_indexed_bam(&[("source-1-90", "chr1", 90), ("source-1-110", "chr1", 110)]);
        let source_2 = create_indexed_bam(&[("source-2-105", "chr1", 105)]);
        let mut reader = multi_reader(&[
            ("source-0.bam", &source_0),
            ("source-1.bam", &source_1),
            ("source-2.bam", &source_2),
        ]);

        reader.fetch(("chr1", 0, 200)).unwrap();

        assert_eq!(
            collect_qnames(&mut reader),
            vec![
                "source-1-90",
                "source-0-100",
                "source-2-105",
                "source-1-110",
                "source-0-120",
            ]
        );
    }

    #[test]
    fn next_source_index_fills_missing_records_before_selecting_minimum() {
        let source_0 = create_indexed_bam(&[("source-0-120", "chr1", 120)]);
        let source_1 = create_indexed_bam(&[("source-1-90", "chr1", 90)]);
        let mut reader = MultiBamReader::new(vec![
            source_reader("source-0.bam", &source_0),
            source_reader("source-1.bam", &source_1),
        ]);
        reader.fetch(("chr1", 0, 200)).unwrap();

        let source_index = reader
            .next_source_index()
            .unwrap()
            .expect("a source should be selected");

        assert_eq!(source_index, 1);
        assert_eq!(reader.sources[0].current_record_sort_key, Some((0, 120, 0)));
        assert_eq!(reader.sources[1].current_record_sort_key, Some((0, 90, 0)));
    }

    #[test]
    fn buffered_source_reader_caches_record_key_and_uses_caller_source_index() {
        let bam_file = create_indexed_bam(&[("source-0-120", "chr1", 120)]);
        let mut source = BufferedSourceReader::new(source_reader("source-0.bam", &bam_file));
        source.reader.fetch(("chr1", 0, 200)).unwrap();

        source.fill_current_record().unwrap();

        assert_eq!(source.current_sort_key(7), Some((0, 120, 7, 0)));
        assert_eq!(source.current_sort_key(3), Some((0, 120, 3, 0)));
    }

    #[test]
    fn multi_source_equal_coordinate_records_are_deterministic() {
        let source_0 = create_indexed_bam(&[
            ("source-0-first", "chr1", 100),
            ("source-0-second", "chr1", 100),
        ]);
        let source_1 = create_indexed_bam(&[
            ("source-1-first", "chr1", 100),
            ("source-1-second", "chr1", 100),
        ]);
        let source_2 = create_indexed_bam(&[("source-2-first", "chr1", 100)]);
        let mut reader = multi_reader(&[
            ("source-0.bam", &source_0),
            ("source-1.bam", &source_1),
            ("source-2.bam", &source_2),
        ]);

        reader.fetch(("chr1", 0, 200)).unwrap();

        assert_eq!(
            collect_qnames(&mut reader),
            vec![
                "source-0-first",
                "source-0-second",
                "source-1-first",
                "source-1-second",
                "source-2-first",
            ]
        );
    }

    #[test]
    fn multi_source_fetch_errors_report_source_and_queried_region() {
        let primary = create_indexed_bam(&[("primary-in-region", "chr1", 100)]);
        let failed = create_indexed_bam(&[("failed-in-region", "chr1", 110)]);
        let mut reader = multi_reader(&[("primary.bam", &primary), ("failed.bam", &failed)]);

        let error = reader.fetch(("missing", 50, 200)).unwrap_err();

        assert!(error.contains("primary.bam"));
        assert!(error.contains("missing:50-200"));
    }

    #[test]
    fn multi_source_records_do_not_receive_source_metadata() {
        let primary = create_indexed_bam(&[("primary-in-region", "chr1", 100)]);
        let mut direct = bam::IndexedReader::from_path(primary.path()).unwrap();
        direct.fetch(("chr1", 50, 200)).unwrap();
        let mut expected = Record::new();
        direct.read(&mut expected).unwrap().unwrap();

        let mut reader = multi_reader(&[("primary.bam", &primary)]);
        reader.fetch(("chr1", 50, 200)).unwrap();
        let mut actual = Record::new();
        reader.read(&mut actual).unwrap().unwrap();

        assert_eq!(actual, expected);
    }
}
