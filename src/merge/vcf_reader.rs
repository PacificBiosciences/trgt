use crate::{merge::tpool::HtsThreadPool, utils::Result};
use clap::ValueEnum;
use rust_htslib::{
    bcf::{self, header::HeaderView, Header, HeaderRecord, Read},
    bgzf,
};
use semver::Version;
use std::{
    collections::{HashMap, HashSet},
    ffi::OsString,
    fs::File,
    io::Read as ReadIo,
    path::{Path, PathBuf},
    sync::Arc,
};

pub const MIN_VERSION_WITH_PADDING_BASE: u64 = 1;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum VcfReadMode {
    Indexed,
    Streaming,
}

impl From<bool> for VcfReadMode {
    fn from(value: bool) -> Self {
        if value {
            VcfReadMode::Streaming
        } else {
            VcfReadMode::Indexed
        }
    }
}

impl std::fmt::Display for VcfReadMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VcfReadMode::Indexed => write!(f, "indexed"),
            VcfReadMode::Streaming => write!(f, "streaming"),
        }
    }
}

enum ReaderBackend {
    Indexed(bcf::IndexedReader),
    Stream(bcf::Reader),
}

impl ReaderBackend {
    fn empty_record(&self) -> bcf::Record {
        match self {
            ReaderBackend::Indexed(reader) => reader.empty_record(),
            ReaderBackend::Stream(reader) => reader.empty_record(),
        }
    }
}

fn add_extension(path: &Path, ext: &str) -> PathBuf {
    let mut out = path.to_path_buf();
    let new_ext: OsString = match path.extension() {
        Some(old) => {
            let mut s = old.to_os_string();
            s.push(".");
            s.push(ext);
            s
        }
        None => OsString::from(ext),
    };
    out.set_extension(new_ext);
    out
}

fn is_indexed_local(file: &Path) -> bool {
    add_extension(file, "csi").exists() || add_extension(file, "tbi").exists()
}

fn has_gzip_magic(path: &Path) -> Result<bool> {
    let mut f = File::open(path).map_err(|e| format!("Failed to open {}: {e}", path.display()))?;
    let mut m = [0u8; 2];
    let n = f
        .read(&mut m)
        .map_err(|e| format!("Failed to read {}: {e}", path.display()))?;
    Ok(n == 2 && m == [0x1f, 0x8b])
}

fn looks_like_bcf(p: &[u8]) -> bool {
    p.len() >= 5 && &p[..3] == b"BCF" && p[3] == 0x02
}

fn looks_like_vcf(p: &[u8]) -> bool {
    let mut p = p;
    // Strip the UTF-8 BOM if present
    if p.starts_with(&[0xEF, 0xBB, 0xBF]) {
        p = &p[3..];
    }
    if p.starts_with(b"##fileformat=VCF") {
        return true;
    }
    false
}

pub fn validate_bgzip_vcf(file: &Path) -> Result<()> {
    let is_bgzf = rust_htslib::bgzf::is_bgzip(file).map_err(|e| {
        format!(
            "Failed to determine whether {} is BGZF-compressed: {e}",
            file.display()
        )
    })?;

    if !is_bgzf {
        if has_gzip_magic(file)? {
            return Err(format!("File {} is gzip-compressed but not BGZF (bgzip). Recompress with bgzip (or `bcftools view -Oz`) and index with `tabix -p vcf` / `bcftools index` as needed.", file.display()));
        }
        return Err(format!(
            "File {} is not BGZF (bgzip) compressed",
            file.display()
        ));
    }

    let mut r = bgzf::Reader::from_path(file)
        .map_err(|e| format!("Failed to open BGZF reader for {}: {e}", file.display()))?;
    let mut buf = vec![0u8; 512];

    let n = r
        .read(&mut buf)
        .map_err(|e| format!("Failed to read BGZF stream from {}: {e}", file.display()))?;
    buf.truncate(n);

    if buf.is_empty() {
        return Err(format!("File {} is empty", file.display()));
    }

    if looks_like_vcf(&buf) || looks_like_bcf(&buf) {
        Ok(())
    } else {
        Err(format!(
            "File {} is BGZF-compressed, but the decompressed header does not look like VCF or BCF",
            file.display()
        ))
    }
}

fn validate_indexed_vcf(file: &Path) -> Result<()> {
    if !is_indexed_local(file) {
        return Err(format!(
            "VCF file {} is not indexed (.tbi or .csi not found)",
            file.display()
        ));
    }
    validate_bgzip_vcf(file)
}

pub struct VcfReader {
    // backend must drop BEFORE _thread_pool to ensure pool is alive during close
    backend: ReaderBackend,
    _thread_pool: Option<Arc<HtsThreadPool>>,
    mode: VcfReadMode,
    pub header: HeaderView,
    current_record: bcf::Record,
    spare_record: Option<bcf::Record>,
    pub needs_pos_adjustment: bool,
    pub index: usize,
    pub sample_n: u32,
    pub file_path: String,
}

impl VcfReader {
    pub fn new(file: PathBuf, index: usize, mode: VcfReadMode) -> Result<Self> {
        Self::new_with_thread_pool(file, index, mode, None)
    }

    pub fn new_with_thread_pool(
        file: PathBuf,
        index: usize,
        mode: VcfReadMode,
        thread_pool: Option<Arc<HtsThreadPool>>,
    ) -> Result<Self> {
        log::debug!("Start opening VCF {:?}", &file.display());

        let (backend, header) = match mode {
            VcfReadMode::Indexed => {
                validate_indexed_vcf(&file)?;
                let mut reader = bcf::IndexedReader::from_path(&file)
                    .map_err(|e| format!("Failed to open VCF file {}: {}", file.display(), e))?;

                if let Some(ref pool) = thread_pool {
                    // SAFETY: Reader is opened and owned here, pool is kept alive alongside
                    // the reader and `attach_to_indexed_reader` only forwards the inner htsFile
                    unsafe {
                        pool.attach_to_indexed_reader(&mut reader)?;
                    }
                }

                let header = reader.header().clone();
                (ReaderBackend::Indexed(reader), header)
            }
            VcfReadMode::Streaming => {
                validate_bgzip_vcf(&file)?;

                let mut reader = bcf::Reader::from_path(&file)
                    .map_err(|e| format!("Failed to open VCF file {}: {}", file.display(), e))?;

                if let Some(ref pool) = thread_pool {
                    // SAFETY: Reader is opened and owned here, pool is kept alive alongside
                    // the reader and `attach_to_reader` only forwards the inner htsFile
                    unsafe {
                        pool.attach_to_reader(&mut reader)?;
                    }
                }

                let header = reader.header().clone();
                (ReaderBackend::Stream(reader), header)
            }
        };

        let version = get_trgt_version(&header, &file)?;
        let sample_n = header.sample_count();
        let needs_pos_adjustment = version.major < MIN_VERSION_WITH_PADDING_BASE;

        log::debug!(
            "{:?} has TRGT version: {}, samples n = {}",
            file.file_name().unwrap(),
            version,
            sample_n
        );

        let current_record = backend.empty_record();

        log::debug!("Finished opening VCF {:?}", &file.display());
        Ok(VcfReader {
            backend,
            _thread_pool: thread_pool,
            header,
            current_record,
            spare_record: None,
            needs_pos_adjustment,
            index,
            sample_n,
            file_path: file.to_string_lossy().into_owned(),
            mode,
        })
    }

    pub fn mode(&self) -> VcfReadMode {
        self.mode
    }

    pub(crate) fn indexed_reader_mut(&mut self) -> Option<&mut bcf::IndexedReader> {
        match &mut self.backend {
            ReaderBackend::Indexed(reader) => Some(reader),
            ReaderBackend::Stream(_) => None,
        }
    }

    pub fn advance(&mut self) -> bool {
        let has_record = match &mut self.backend {
            ReaderBackend::Indexed(reader) => {
                Self::advance_reader(reader, &mut self.current_record)
            }
            ReaderBackend::Stream(reader) => Self::advance_reader(reader, &mut self.current_record),
        };
        if has_record {
            self.update_record_for_version();
        }
        has_record
    }

    fn advance_reader<R: Read>(reader: &mut R, record: &mut bcf::Record) -> bool {
        match reader.read(record) {
            Some(Ok(())) => true,
            Some(Err(_)) | None => false,
        }
    }

    pub fn take_current_record(&mut self) -> bcf::Record {
        let mut next = self
            .spare_record
            .take()
            .unwrap_or_else(|| self.backend.empty_record());
        next.clear();
        std::mem::swap(&mut self.current_record, &mut next);
        next
    }

    pub fn return_record(&mut self, record: bcf::Record) {
        record.clear();
        self.spare_record = Some(record);
    }

    fn update_record_for_version(&mut self) {
        if !self.needs_pos_adjustment {
            return;
        }

        let al_fmt = self
            .current_record
            .format(b"AL")
            .integer()
            .expect("Error accessing FORMAT AL");
        let al_0 = *al_fmt[0].iter().min().unwrap();
        if al_0 != 0 {
            self.current_record.set_pos(self.current_record.pos() - 1);
        }
    }
}

pub fn get_trgt_version(vcf_header: &HeaderView, file: &Path) -> Result<Version> {
    let mut trgt_version = None;

    for record in vcf_header.header_records().iter() {
        if let HeaderRecord::Generic { key, value } = record {
            if key == "trgtVersion" {
                trgt_version = Some(value.clone());
            }
        }
    }

    if trgt_version.is_none() {
        let mut has_allr = false;
        let mut has_alci = false;
        let mut is_integer_am = false;
        for record in vcf_header.header_records().iter() {
            if let HeaderRecord::Format { key: _, values } = record {
                if let Some(id) = values.get("ID") {
                    match id.as_str() {
                        "ALLR" => has_allr = true,
                        "ALCI" => has_alci = true,
                        "AM" => {
                            if let Some(typ) = values.get("Type") {
                                is_integer_am = typ == "Integer";
                            }
                        }
                        _ => {}
                    }
                }
            }
        }

        if has_alci {
            trgt_version = Some("0.3.4".to_string());
        } else if has_allr && is_integer_am {
            trgt_version = Some("0.4.0".to_string());
        }

        if trgt_version.is_none() {
            return Err(format!("Non-TRGT VCF supplied {}", file.display()));
        }
    }

    let version = Version::parse(&trgt_version.unwrap())
        .map_err(|e| format!("Failed to parse version: {}", e))?;

    Ok(version)
}

pub struct VcfReaders {
    // readers must drop BEFORE _thread_pool to ensure pool is alive during close
    pub readers: Vec<VcfReader>,
    _thread_pool: Option<Arc<HtsThreadPool>>,
    mode: VcfReadMode,
}

impl VcfReaders {
    pub fn new(
        vcf_files: Vec<PathBuf>,
        mode: VcfReadMode,
        thread_pool: Option<Arc<HtsThreadPool>>,
    ) -> Result<Self> {
        let readers = vcf_files
            .into_iter()
            .enumerate()
            .map(|(index, file)| {
                VcfReader::new_with_thread_pool(file, index, mode, thread_pool.clone())
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(VcfReaders {
            readers,
            _thread_pool: thread_pool,
            mode,
        })
    }

    pub fn mode(&self) -> VcfReadMode {
        self.mode
    }

    pub fn get_contig_order(&self) -> Result<Vec<String>> {
        match self.mode {
            VcfReadMode::Indexed => self.collect_contig_union(),
            VcfReadMode::Streaming => self.verify_identical_contig_order(),
        }
    }

    fn collect_contig_union(&self) -> Result<Vec<String>> {
        let mut contig_map: HashMap<String, HashSet<u64>> = HashMap::new();
        let mut contig_order = Vec::new();
        for reader in &self.readers {
            for record in reader.header.header_records() {
                if let HeaderRecord::Contig { values, .. } = record {
                    let id = values.get("ID").unwrap().to_string();
                    let length = values
                        .get("length")
                        .and_then(|l| l.parse::<u64>().ok())
                        .unwrap_or(0);
                    let entry = contig_map.entry(id.clone()).or_insert_with(|| {
                        contig_order.push(id.clone());
                        HashSet::new()
                    });
                    entry.insert(length);
                }
            }
        }
        for id in &contig_order {
            let lengths = contig_map.get(id).unwrap();
            if lengths.len() > 1 {
                return Err(format!(
                    "Inconsistent contig definitions found in VCF headers: Contig '{}' is defined with multiple lengths: {:?}",
                    id, lengths
                ));
            }
        }
        Ok(contig_order)
    }

    fn verify_identical_contig_order(&self) -> Result<Vec<String>> {
        let mut expected_contigs: Option<Vec<(String, Option<u64>)>> = None;

        for reader in &self.readers {
            let contigs = contigs_in_header_order(&reader.header, &reader.file_path)?;
            if let Some(ref expected) = expected_contigs {
                if contigs.len() != expected.len() {
                    return Err(format!(
                        "Contig count mismatch for {}: expected {} contigs, found {}",
                        reader.file_path,
                        expected.len(),
                        contigs.len()
                    ));
                }

                for (idx, ((id, length), (expected_id, expected_length))) in
                    contigs.iter().zip(expected.iter()).enumerate()
                {
                    if id != expected_id {
                        return Err(format!(
                            "Contig order differs at position {}: expected '{}' but found '{}' in {}",
                            idx + 1,
                            expected_id,
                            id,
                            reader.file_path
                        ));
                    }

                    if let (Some(current_len), Some(expected_len)) = (length, expected_length) {
                        if current_len != expected_len {
                            return Err(format!(
                                "Contig length mismatch for '{}' in {}: expected {}, found {}",
                                id, reader.file_path, expected_len, current_len
                            ));
                        }
                    }
                }
            } else {
                expected_contigs = Some(contigs);
            }
        }

        Ok(expected_contigs
            .unwrap_or_default()
            .into_iter()
            .map(|(id, _)| id)
            .collect())
    }

    pub fn merge_headers(&self, dst_header: &mut Header) -> Result<()> {
        let mut observed_sample_ids = HashSet::new();

        for reader in &self.readers {
            let src_header = &reader.header;
            // SAFETY: All headers come from rust-htslib and stay alive during merge, the merge
            // function only mutates the destination header
            unsafe {
                dst_header.inner =
                    rust_htslib::htslib::bcf_hdr_merge(dst_header.inner, src_header.inner);
            }

            for sample_id in src_header.samples() {
                if observed_sample_ids.contains(sample_id) {
                    return Err(format!(
                        "Duplicate sample ID found: {}",
                        String::from_utf8_lossy(sample_id)
                    ));
                }
                observed_sample_ids.insert(sample_id.to_vec());
                dst_header.push_sample(sample_id);
            }
        }

        // SAFETY: `dst_header.inner` remains valid; sync refreshes derived header state after the merge
        unsafe {
            rust_htslib::htslib::bcf_hdr_sync(dst_header.inner);
        }

        Ok(())
    }

    pub fn len(&self) -> usize {
        self.readers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.readers.is_empty()
    }

    #[cfg(test)]
    pub fn empty(mode: VcfReadMode) -> Self {
        VcfReaders {
            readers: Vec::new(),
            _thread_pool: None,
            mode,
        }
    }
}

fn contigs_in_header_order(
    header: &HeaderView,
    file_path: &str,
) -> Result<Vec<(String, Option<u64>)>> {
    let mut contigs = Vec::new();
    for record in header.header_records() {
        if let HeaderRecord::Contig { values, .. } = record {
            let id = values
                .get("ID")
                .ok_or_else(|| format!("Contig header missing ID in {}", file_path))?;
            let length = values.get("length").and_then(|l| l.parse::<u64>().ok());
            contigs.push((id.to_string(), length));
        }
    }
    Ok(contigs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_util::{TestVcfBuilder, TestVcfRecord};
    use rust_htslib::bcf::{index, record::GenotypeAllele};
    use tempfile::tempdir;

    #[test]
    fn test_add_extension() {
        let path = PathBuf::from("file.vcf");
        let ext = "gz";
        let new_path = add_extension(&path, ext);
        assert_eq!(new_path, PathBuf::from("file.vcf.gz"));

        let path_with_dir = PathBuf::from("/a/b/file.vcf");
        let new_path_with_dir = add_extension(&path_with_dir, ext);
        assert_eq!(new_path_with_dir, PathBuf::from("/a/b/file.vcf.gz"));
    }

    #[test]
    fn streaming_reader_accepts_unindexed_vcf() -> Result<()> {
        let temp_vcf = TestVcfBuilder::new()
            .header_line("trgtVersion=1.0.0")
            .contig("chr1", 100)
            .sample("sample1")
            .record(
                TestVcfRecord::new()
                    .pos(10)
                    .alleles(&["A", "AT"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build_with_compression(true);

        let readers = VcfReaders::new(
            vec![temp_vcf.path().to_path_buf()],
            VcfReadMode::Streaming,
            None,
        )?;
        assert_eq!(readers.len(), 1);
        assert_eq!(readers.mode(), VcfReadMode::Streaming);
        Ok(())
    }

    #[test]
    fn indexed_reader_requires_index() {
        let temp_vcf = TestVcfBuilder::new()
            .header_line("trgtVersion=1.0.0")
            .contig("chr1", 100)
            .sample("sample1")
            .record(
                TestVcfRecord::new()
                    .pos(10)
                    .alleles(&["A", "AT"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build_with_compression(true);

        let result = VcfReaders::new(
            vec![temp_vcf.path().to_path_buf()],
            VcfReadMode::Indexed,
            None,
        );
        assert!(result.is_err());
    }

    #[test]
    fn streaming_contig_order_mismatch_errors() -> Result<()> {
        let first = TestVcfBuilder::new()
            .header_line("trgtVersion=1.0.0")
            .contig("chr1", 100)
            .contig("chr2", 200)
            .sample("s1")
            .record(
                TestVcfRecord::new()
                    .pos(10)
                    .alleles(&["A", "AT"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build_with_compression(true);

        let second = TestVcfBuilder::new()
            .header_line("trgtVersion=1.0.0")
            .contig("chr2", 200)
            .contig("chr1", 100)
            .sample("s1")
            .record(
                TestVcfRecord::new()
                    .pos(10)
                    .alleles(&["A", "AT"])
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
            )
            .build_with_compression(true);

        let readers = VcfReaders::new(
            vec![first.path().to_path_buf(), second.path().to_path_buf()],
            VcfReadMode::Streaming,
            None,
        )?;

        let contigs = readers.get_contig_order();
        assert!(contigs.is_err());
        Ok(())
    }

    #[test]
    fn indexed_contig_order_allows_union() -> Result<()> {
        let dir = tempdir().unwrap();
        let path_a = dir.path().join("a.vcf.gz");
        let path_b = dir.path().join("b.vcf.gz");

        let build_indexed = |path: &Path| {
            let file = TestVcfBuilder::new()
                .header_line("trgtVersion=1.0.0")
                .contig("chr1", 100)
                .contig("chr2", 200)
                .sample("s1")
                .record(
                    TestVcfRecord::new()
                        .pos(10)
                        .alleles(&["A", "AT"])
                        .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]),
                )
                .build_with_compression(true);
            std::fs::copy(file.path(), path).unwrap();
            index::build(path, None::<&Path>, 1, index::Type::Tbx).map_err(|e| e.msg)
        };

        build_indexed(&path_a)?;
        build_indexed(&path_b)?;

        let readers = VcfReaders::new(vec![path_a, path_b], VcfReadMode::Indexed, None)?;

        let contigs = readers.get_contig_order()?;
        assert_eq!(contigs, vec!["chr1".to_string(), "chr2".to_string()]);
        Ok(())
    }
}
