use rust_htslib::{
    bam::{
        self,
        header::HeaderRecord as BamHeaderRecord,
        record::{Aux, Cigar, CigarString},
        Record as BamRecord, Writer as BamWriter,
    },
    bcf::{header::Header, record::GenotypeAllele, Format, Writer},
};
use std::{
    fs,
    path::{Path, PathBuf},
};
use tempfile::{NamedTempFile, TempDir};

use super::InputSource;

#[derive(Default)]
pub struct TestVcfBuilder {
    contigs: Vec<(String, u64)>,
    header_lines: Vec<String>,
    samples: Vec<String>,
    records: Vec<TestVcfRecord>,
}

#[derive(Default)]
pub struct TestVcfRecord {
    rid: Option<u32>,
    pos: i64,
    id: Option<Vec<u8>>,
    qual: Option<f32>,
    filters: Vec<String>,
    alleles: Vec<Vec<u8>>,
    info: Vec<(String, InfoValue)>,
    genotypes: Vec<GenotypeAllele>,
    format: Vec<(String, FormatValue)>,
}

enum InfoValue {
    String(Vec<u8>),
    Integer(Vec<i32>),
    Flag,
}

enum FormatValue {
    String(Vec<Vec<u8>>),
    Integer(Vec<i32>),
    Float(Vec<f32>),
}

pub struct TestBam {
    _dir: TempDir,
    path: PathBuf,
}

impl TestBam {
    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn source(&self) -> InputSource {
        InputSource::Local(self.path.clone())
    }
}

pub fn create_indexed_bam(records: &[(&str, &str, i64)]) -> NamedTempFile {
    create_indexed_bam_from_records(
        records
            .iter()
            .map(|(name, contig, pos)| ((*name).to_owned(), (*contig).to_owned(), *pos)),
    )
}

pub fn create_indexed_bam_from_records<I>(records: I) -> NamedTempFile
where
    I: IntoIterator<Item = (String, String, i64)>,
{
    let bam_file = tempfile::Builder::new()
        .suffix(".bam")
        .tempfile()
        .expect("failed to create temporary BAM path");

    let mut header = bam::Header::new();
    for contig in ["chr1", "chr2"] {
        let mut sq = BamHeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", 1_000);
        header.push_record(&sq);
    }

    {
        let mut writer = BamWriter::from_path(bam_file.path(), &header, bam::Format::Bam)
            .expect("failed to create temporary BAM");
        for (name, contig, pos) in records {
            writer
                .write(&mapped_bam_record(&name, test_bam_contig_tid(&contig), pos))
                .expect("failed to write temporary BAM record");
        }
    }

    bam::index::build(bam_file.path(), None, bam::index::Type::Bai, 1)
        .expect("failed to index temporary BAM");

    bam_file
}

pub fn create_indexed_test_bam(name: &str, contigs: &[(&str, u32)], samples: &[&str]) -> TestBam {
    let bam = create_test_bam(name, contigs, samples);
    bam::index::build(bam.path(), None, bam::index::Type::Bai, 1)
        .expect("failed to index temporary BAM");
    bam
}

pub fn create_unindexed_test_bam(name: &str, contigs: &[(&str, u32)], samples: &[&str]) -> TestBam {
    create_test_bam(name, contigs, samples)
}

fn create_test_bam(name: &str, contigs: &[(&str, u32)], samples: &[&str]) -> TestBam {
    let temp_dir = tempfile::tempdir().expect("failed to create temp dir");
    let path = temp_dir.path().join(name);
    write_test_bam(&path, contigs, samples);
    TestBam {
        _dir: temp_dir,
        path,
    }
}

fn write_test_bam(path: &Path, contigs: &[(&str, u32)], samples: &[&str]) {
    let mut header = bam::Header::new();
    for &(contig, len) in contigs {
        let mut sq = BamHeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", len);
        header.push_record(&sq);
    }
    for (index, &sample) in samples.iter().enumerate() {
        let mut rg = BamHeaderRecord::new(b"RG");
        rg.push_tag(b"ID", format!("rg{index}"));
        rg.push_tag(b"SM", sample);
        header.push_record(&rg);
    }

    let mut writer = BamWriter::from_path(path, &header, bam::Format::Bam)
        .expect("failed to create temporary BAM");
    if !contigs.is_empty() {
        writer
            .write(&mapped_bam_record("read1", 0, 10))
            .expect("failed to write temporary BAM record");
    }
}

fn test_bam_contig_tid(contig: &str) -> i32 {
    match contig {
        "chr1" => 0,
        "chr2" => 1,
        _ => panic!("unexpected contig {contig}"),
    }
}

fn mapped_bam_record(name: &str, tid: i32, pos: i64) -> BamRecord {
    let bases = b"ACGTACGT";
    let quals = vec![40; bases.len()];
    let mut record = BamRecord::new();
    record.set(
        name.as_bytes(),
        Some(&CigarString(vec![Cigar::Match(bases.len() as u32)])),
        bases,
        &quals,
    );
    record.set_tid(tid);
    record.set_pos(pos);
    record.set_mapq(60);
    record.unset_unmapped();
    record.push_aux(b"rq", Aux::Float(0.99)).unwrap();
    record
}

impl TestVcfRecord {
    pub fn new() -> Self {
        Self {
            rid: Some(0),
            alleles: vec![b"A".to_vec()],
            ..Default::default()
        }
    }

    pub fn pos(mut self, pos: i64) -> Self {
        self.pos = pos;
        self
    }

    pub fn id<T: AsRef<[u8]>>(mut self, id: T) -> Self {
        self.id = Some(id.as_ref().to_vec());
        self
    }

    pub fn qual(mut self, qual: f32) -> Self {
        self.qual = Some(qual);
        self
    }

    pub fn filter<S: ToString>(mut self, filter: S) -> Self {
        self.filters.push(filter.to_string());
        self
    }

    pub fn alleles<T: AsRef<[u8]>>(mut self, alleles: &[T]) -> Self {
        self.alleles = alleles.iter().map(|a| a.as_ref().to_vec()).collect();
        self
    }

    pub fn genotype(mut self, genotypes: &[GenotypeAllele]) -> Self {
        self.genotypes = genotypes.to_vec();
        self
    }

    pub fn info_string<K: ToString, V: AsRef<[u8]>>(mut self, key: K, value: V) -> Self {
        self.info
            .push((key.to_string(), InfoValue::String(value.as_ref().to_vec())));
        self
    }

    pub fn info_integer<K: ToString>(mut self, key: K, values: &[i32]) -> Self {
        self.info
            .push((key.to_string(), InfoValue::Integer(values.to_vec())));
        self
    }

    pub fn info_flag<K: ToString>(mut self, key: K) -> Self {
        self.info.push((key.to_string(), InfoValue::Flag));
        self
    }

    pub fn format_integer<K: ToString>(mut self, key: K, values: &[i32]) -> Self {
        self.format
            .push((key.to_string(), FormatValue::Integer(values.to_vec())));
        self
    }

    pub fn format_float<K: ToString>(mut self, key: K, values: &[f32]) -> Self {
        self.format
            .push((key.to_string(), FormatValue::Float(values.to_vec())));
        self
    }

    pub fn format_string<K: ToString, V: AsRef<[u8]>>(mut self, key: K, values: &[V]) -> Self {
        let converted = values.iter().map(|v| v.as_ref().to_vec()).collect();
        self.format
            .push((key.to_string(), FormatValue::String(converted)));
        self
    }

    pub fn with_trgt_info<S: AsRef<[u8]>>(trid: S, end: i32, motifs: S, struc: S) -> Self {
        Self::new()
            .info_string("TRID", trid)
            .info_integer("END", &[end])
            .info_string("MOTIFS", motifs)
            .info_string("STRUC", struc)
    }

    pub fn with_trgt_format_values(self) -> Self {
        self.genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .format_integer("AL", &[1, 2])
            .format_string("ALLR", &["1-1,1-2"])
            .format_integer("SD", &[5, 6])
            .format_string("MC", &["1,2"])
            .format_string("MS", &["0-1,0-2"])
            .format_float("AP", &[0.9, 0.8])
            .format_float("AM", &[0.1, 0.2])
    }

    pub fn default_trgt<S: AsRef<[u8]>>(trid: S, pos: i64, end: i32) -> Self {
        Self::with_trgt_info(trid.as_ref(), end, b"A".as_ref(), b"(A)n".as_ref())
            .pos(pos)
            .alleles(&["A", "AA"])
            .with_trgt_format_values()
    }
}

impl TestVcfBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn contig<S: ToString>(mut self, name: S, length: u64) -> Self {
        self.contigs.push((name.to_string(), length));
        self
    }

    pub fn sample<S: ToString>(mut self, name: S) -> Self {
        self.samples.push(name.to_string());
        self
    }

    pub fn record(mut self, record: TestVcfRecord) -> Self {
        self.records.push(record);
        self
    }

    pub fn header_line<S: ToString>(mut self, line: S) -> Self {
        self.header_lines.push(line.to_string());
        self
    }

    pub fn add_info<S: ToString>(mut self, id: S, num: &str, type_: &str, desc: &str) -> Self {
        self.header_lines.push(format!(
            "INFO=<ID={},Number={},Type={},Description=\"{}\">",
            id.to_string(),
            num,
            type_,
            desc
        ));
        self
    }

    pub fn add_format<S: ToString>(mut self, id: S, num: &str, type_: &str, desc: &str) -> Self {
        self.header_lines.push(format!(
            "FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
            id.to_string(),
            num,
            type_,
            desc
        ));
        self
    }

    pub fn add_filter<S: ToString>(mut self, id: S, desc: &str) -> Self {
        self.header_lines.push(format!(
            "FILTER=<ID={},Description=\"{}\">",
            id.to_string(),
            desc
        ));
        self
    }

    pub fn build_header(&self) -> Header {
        let mut header = Header::new();
        header.push_record(br#"##fileformat=VCFv4.3"#);

        for (name, length) in &self.contigs {
            header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
        }

        for line in &self.header_lines {
            header.push_record(format!("##{}", line).as_bytes());
        }

        header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);

        for sample in &self.samples {
            header.push_sample(sample.as_bytes());
        }

        header
    }

    pub fn build(self) -> NamedTempFile {
        self.build_with_compression(true)
    }

    pub fn build_with_compression(self, compress: bool) -> NamedTempFile {
        let header = self.build_header();
        let temp_file = NamedTempFile::new().expect("Failed to create temp file");
        let uncompressed = !compress;
        let mut writer = Writer::from_path(temp_file.path(), &header, uncompressed, Format::Vcf)
            .expect("Failed to create writer");

        for rec in &self.records {
            let mut record = writer.empty_record();

            if let Some(rid) = rec.rid {
                record.set_rid(Some(rid));
            }
            record.set_pos(rec.pos);
            if let Some(q) = rec.qual {
                record.set_qual(q);
            }
            if let Some(ref id) = rec.id {
                record.set_id(id).unwrap();
            }

            let allele_refs: Vec<&[u8]> = rec.alleles.iter().map(|a| a.as_slice()).collect();
            record.set_alleles(&allele_refs).unwrap();

            for filter in &rec.filters {
                record.push_filter(filter.as_bytes()).unwrap();
            }

            for (key, value) in &rec.info {
                let key_bytes = key.as_bytes();
                match value {
                    InfoValue::String(v) => record.push_info_string(key_bytes, &[v]).unwrap(),
                    InfoValue::Integer(v) => record.push_info_integer(key_bytes, v).unwrap(),
                    InfoValue::Flag => record.push_info_flag(key_bytes).unwrap(),
                }
            }

            if !rec.genotypes.is_empty() {
                record.push_genotypes(&rec.genotypes).unwrap();
            }

            for (key, value) in &rec.format {
                let key_bytes = key.as_bytes();
                match value {
                    FormatValue::String(v) => {
                        let refs: Vec<&[u8]> = v.iter().map(|x| x.as_slice()).collect();
                        record.push_format_string(key_bytes, &refs).unwrap();
                    }
                    FormatValue::Integer(v) => record.push_format_integer(key_bytes, v).unwrap(),
                    FormatValue::Float(v) => record.push_format_float(key_bytes, v).unwrap(),
                }
            }

            writer.write(&record).unwrap();
        }

        temp_file
    }

    pub fn with_trgt_defaults(self) -> Self {
        self.add_info("TRID", "1", "String", "TR ID")
            .add_info("END", "1", "Integer", "End position")
            .add_info("MOTIFS", "1", "String", "Motifs")
            .add_info("STRUC", "1", "String", "Structure")
    }

    pub fn with_trgt_format_defaults(self) -> Self {
        self.add_format("AL", ".", "Integer", "Allele length")
            .add_format("ALLR", ".", "String", "Allele range")
            .add_format("SD", ".", "Integer", "Spanning depth")
            .add_format("MC", ".", "String", "Motif counts")
            .add_format("MS", ".", "String", "Motif spans")
            .add_format("AP", ".", "Float", "Allele purity")
            .add_format("AM", ".", "Float", "Allele methylation")
    }

    pub fn with_trgt_merge_defaults(self) -> Self {
        self.header_line("trgtVersion=1.0.0")
            .with_trgt_defaults()
            .with_trgt_format_defaults()
    }
}

pub fn write_vcf_with_gt(
    trid: &[u8],
    alleles: &[&[u8]],
    genotype: &[GenotypeAllele],
) -> NamedTempFile {
    TestVcfBuilder::new()
        .contig("chr1", 1)
        .add_info("TRID", "1", "String", "TR id")
        .sample("sample")
        .record(
            TestVcfRecord::new()
                .alleles(alleles)
                .info_string("TRID", trid)
                .genotype(genotype),
        )
        .build()
}

pub fn write_trgt_vcf(
    path: &Path,
    sample: &str,
    trid: &str,
    pos: i64,
) -> std::result::Result<(), String> {
    write_trgt_vcf_with_records(
        vec![TestVcfRecord::default_trgt(trid, pos, (pos + 5) as i32)],
        path,
        sample,
    )
}

pub fn write_trgt_vcf_with_records(
    records: Vec<TestVcfRecord>,
    path: &Path,
    sample: &str,
) -> std::result::Result<(), String> {
    let builder = TestVcfBuilder::new()
        .with_trgt_merge_defaults()
        .contig("chr1", 10000)
        .sample(sample);
    let source = records
        .into_iter()
        .fold(builder, TestVcfBuilder::record)
        .build();
    fs::copy(source.path(), path).map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as BamRead;
    use rust_htslib::bcf::Read;

    #[test]
    fn test_vcf_builder_basic() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("test_sample")
            .record(
                TestVcfRecord::new()
                    .pos(100)
                    .alleles(&["A", "AT"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build();

        assert!(temp_vcf.path().exists());
    }

    #[test]
    fn test_indexed_bam_builder_writes_fetchable_records() {
        let temp_bam = create_indexed_bam(&[
            ("before-region", "chr1", 10),
            ("in-region", "chr1", 100),
            ("other-contig", "chr2", 100),
        ]);

        let mut reader = rust_htslib::bam::IndexedReader::from_path(temp_bam.path()).unwrap();
        reader.fetch(("chr1", 50, 150)).unwrap();

        let mut record = rust_htslib::bam::Record::new();
        let result = reader.read(&mut record).unwrap();

        assert!(result.is_ok());
        assert_eq!(std::str::from_utf8(record.qname()).unwrap(), "in-region");
        assert!(reader.read(&mut record).is_none());
    }

    #[test]
    fn test_indexed_test_bam_builder_writes_sample_metadata() {
        let temp_bam = create_indexed_test_bam("sampled.bam", &[("chr1", 1000)], &["sample-a"]);
        let source = temp_bam.source();
        let header = crate::utils::get_bam_header(&source).unwrap();

        assert_eq!(
            crate::utils::get_sample_name(&source, &header).unwrap(),
            "sample-a"
        );

        let mut reader = rust_htslib::bam::IndexedReader::from_path(temp_bam.path()).unwrap();
        reader.fetch(("chr1", 0, 100)).unwrap();

        let mut record = rust_htslib::bam::Record::new();
        let result = reader.read(&mut record).unwrap();

        assert!(result.is_ok());
        assert_eq!(std::str::from_utf8(record.qname()).unwrap(), "read1");
    }

    #[test]
    fn test_trgt_vcf_helpers_write_default_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("trgt.vcf.gz");

        write_trgt_vcf_with_records(
            vec![TestVcfRecord::default_trgt("repeat1", 1000, 1050)],
            &path,
            "S1",
        )
        .unwrap();

        let mut reader = rust_htslib::bcf::Reader::from_path(path).unwrap();
        assert_eq!(reader.header().sample_count(), 1);

        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let trid = record.info(b"TRID").string().unwrap().unwrap();
        assert_eq!(trid[0], b"repeat1");

        let allr = record.format(b"ALLR").string().unwrap();
        assert_eq!(allr[0], b"1-1,1-2");
    }

    #[test]
    fn test_vcf_builder_with_info() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .add_info("TRID", "1", "String", "TR id")
            .add_info("END", "1", "Integer", "End position")
            .sample("sample")
            .record(
                TestVcfRecord::new()
                    .alleles(&["A", "AA"])
                    .info_string("TRID", "my_trid")
                    .info_integer("END", &[200])
                    .genotype(&[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)]),
            )
            .build();
        assert!(temp_vcf.path().exists());
    }

    #[test]
    fn test_vcf_builder_with_qual() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("sample")
            .record(
                TestVcfRecord::new()
                    .pos(100)
                    .qual(42.0)
                    .alleles(&["A", "T"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();
        assert!((record.qual() - 42.0).abs() < 0.001);
    }

    #[test]
    fn test_vcf_builder_with_format_fields() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .add_format("AL", ".", "Integer", "Allele length")
            .add_format("AP", ".", "Float", "Allele purity")
            .sample("S1")
            .sample("S2")
            .record(
                TestVcfRecord::new()
                    .pos(100)
                    .alleles(&["A", "AT", "ATT"])
                    .genotype(&[
                        GenotypeAllele::Unphased(0),
                        GenotypeAllele::Unphased(1),
                        GenotypeAllele::Unphased(1),
                        GenotypeAllele::Unphased(2),
                    ])
                    .format_integer("AL", &[10, 15, 15, 20])
                    .format_float("AP", &[0.95, 0.98, 0.97, 0.99]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let al = record.format(b"AL").integer().unwrap();
        assert_eq!(al[0], &[10, 15]);
        assert_eq!(al[1], &[15, 20]);

        let ap = record.format(b"AP").float().unwrap();
        assert!((ap[0][0] - 0.95).abs() < 0.001);
        assert!((ap[1][1] - 0.99).abs() < 0.001);
    }

    #[test]
    fn test_vcf_builder_multiple_records() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .add_info("TRID", "1", "String", "TR id")
            .add_info("END", "1", "Integer", "End position")
            .add_info("MOTIFS", "1", "String", "Motifs")
            .add_info("STRUC", "1", "String", "Structure")
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .info_string("TRID", "repeat1")
                    .info_integer("END", &[1020])
                    .info_string("MOTIFS", "T")
                    .info_string("STRUC", "(T)n"),
            )
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .info_string("TRID", "repeat2")
                    .info_integer("END", &[1020])
                    .info_string("MOTIFS", "T")
                    .info_string("STRUC", "(T)n"),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut count = 0;
        let mut trids = Vec::new();
        for result in reader.records() {
            let record = result.unwrap();
            let trid = record.info(b"TRID").string().unwrap().unwrap();
            trids.push(String::from_utf8_lossy(trid[0]).to_string());
            count += 1;
        }
        assert_eq!(count, 2);
        assert_eq!(trids, vec!["repeat1", "repeat2"]);
    }

    #[test]
    fn test_build_header_only() {
        let builder = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .add_info("TRID", "1", "String", "TR id")
            .add_format("AL", ".", "Integer", "Allele length")
            .sample("S1")
            .sample("S2");

        let header = builder.build_header();
        let temp_file = NamedTempFile::new().unwrap();
        let writer = Writer::from_path(temp_file.path(), &header, false, Format::Vcf).unwrap();
        assert_eq!(writer.header().sample_count(), 2);
    }

    #[test]
    fn test_write_vcf_with_gt_backwards_compat() {
        let temp_vcf = write_vcf_with_gt(
            b"test_trid",
            &[b"A", b"AA"],
            &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
        );
        assert!(temp_vcf.path().exists());
    }

    #[test]
    fn test_vcf_builder_with_id() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("sample")
            .record(
                TestVcfRecord::new()
                    .pos(100)
                    .id("rs12345")
                    .alleles(&["A", "T"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();
        let id = record.id();
        assert_eq!(id, b"rs12345");
    }

    #[test]
    fn test_vcf_builder_with_filter() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .add_filter("LowQual", "Low quality")
            .sample("sample")
            .record(
                TestVcfRecord::new()
                    .pos(100)
                    .alleles(&["A", "T"])
                    .filter("LowQual")
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();
        assert!(record.has_filter("LowQual".as_bytes()));
    }

    #[test]
    fn test_vcf_builder_with_info_flag() {
        let temp_vcf = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .add_info("TestFlag", "0", "Flag", "Test flag")
            .sample("sample")
            .record(
                TestVcfRecord::new()
                    .pos(100)
                    .alleles(&["A", "T"])
                    .info_flag("TestFlag")
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();
        let flag = record.info(b"TestFlag").flag().unwrap();
        assert!(flag);
    }

    #[test]
    fn test_vcf_builder_trgt_preset() {
        let temp_vcf = TestVcfBuilder::new()
            .with_trgt_defaults()
            .contig("chr1", 10000)
            .sample("sample")
            .record(
                TestVcfRecord::with_trgt_info("repeat1", 1020, "CAG", "(CAG)n")
                    .pos(999)
                    .alleles(&["A", "ACAG"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build();

        let mut reader = rust_htslib::bcf::Reader::from_path(temp_vcf.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let trid = record.info(b"TRID").string().unwrap().unwrap();
        assert_eq!(trid[0], b"repeat1");

        let end = record.info(b"END").integer().unwrap().unwrap();
        assert_eq!(end[0], 1020);

        let motifs = record.info(b"MOTIFS").string().unwrap().unwrap();
        assert_eq!(motifs[0], b"CAG");

        let struc = record.info(b"STRUC").string().unwrap().unwrap();
        assert_eq!(struc[0], b"(CAG)n");
    }
}
