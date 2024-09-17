use crate::utils::Result;
use rust_htslib::bcf::{self, header::HeaderView, Header, HeaderRecord, Read};
use semver::Version;
use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader, Seek},
    path::{Path, PathBuf},
};

const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

fn add_extension(path: &Path, ext: &str) -> PathBuf {
    let mut file_name = path.file_name().unwrap().to_os_string();
    file_name.push(".");
    file_name.push(ext);
    path.with_file_name(file_name)
}

fn is_indexed(file: &Path) -> bool {
    let csi_index = add_extension(file, "csi");
    let tbi_index = add_extension(file, "tbi");
    csi_index.exists() || tbi_index.exists()
}

fn validate_vcf_file(file: &Path) -> Result<()> {
    let mut fp = std::fs::File::open(file).map_err(|e| format!("Failed to open file: {}", e))?;
    let mut buffer = [0; 2];
    std::io::Read::read_exact(&mut fp, &mut buffer)
        .map_err(|e| format!("Failed to read file: {}", e))?;

    if buffer != GZIP_MAGIC_NUMBER {
        return Err(format!("File {} is not gzip compressed", file.display()));
    }

    fp.rewind()
        .map_err(|e| format!("Failed to rewind file: {}", e))?;

    let gz = flate2::read::GzDecoder::new(fp);
    let mut reader = BufReader::new(gz);
    let mut first_line = String::new();
    reader
        .read_line(&mut first_line)
        .map_err(|e| format!("Failed to read file: {}", e))?;

    if !first_line.trim().starts_with("##fileformat=VCFv") {
        return Err(format!("File {} is not a valid VCF file", file.display()));
    }

    Ok(())
}

pub struct VcfReader {
    pub reader: bcf::IndexedReader,
    pub header: bcf::header::HeaderView,
    pub current_record: bcf::Record,
    pub version: Version,
    pub index: usize,
    pub sample_n: u32,
    pub file_path: String,
}

impl VcfReader {
    pub fn new(file: PathBuf, index: usize) -> Result<Self> {
        log::debug!("Start opening VCF {:?}", &file);

        validate_vcf_file(&file).map_err(|e| format!("Error validating VCF: {}", e))?;

        if !is_indexed(&file) {
            return Err(format!("VCF file {} is not indexed", file.display()));
        }

        let reader = bcf::IndexedReader::from_path(&file)
            .map_err(|e| format!("Failed to open VCF file {}: {}", file.display(), e))?;
        let header = reader.header().clone();

        let version = get_trgt_version(&header, &file)?;
        let sample_n = header.sample_count();

        log::debug!(
            "{:?} has version: {}, samples n = {}",
            file.file_name().unwrap(),
            version,
            sample_n
        );

        let current_record = reader.empty_record();
        log::debug!("Finished opening VCF {:?}", &file);
        Ok(VcfReader {
            reader,
            header,
            current_record,
            version,
            index,
            sample_n,
            file_path: file.to_string_lossy().into_owned(),
        })
    }

    pub fn advance(&mut self) -> bool {
        match self.reader.read(&mut self.current_record) {
            Some(Ok(())) => {
                self.update_record_for_version();
                true
            }
            Some(Err(_)) | None => false,
        }
    }
    fn update_record_for_version(&mut self) {
        if self.version.major < Version::new(1, 0, 0).major {
            // Only zero-length alleles had padding in earlier versions
            let al_0 = self
                .current_record
                .format(b"AL")
                .integer()
                .expect("Error accessing FORMAT AL")[0]
                .iter()
                .min()
                .cloned()
                .unwrap();
            if al_0 != 0 {
                self.current_record.set_pos(self.current_record.pos() - 1);
            }
        }
    }
}

fn get_trgt_version(vcf_header: &HeaderView, file: &Path) -> Result<Version> {
    let mut trgt_version = None;

    for record in vcf_header.header_records().iter() {
        if let HeaderRecord::Generic { key, value } = record {
            if key == "trgtVersion" {
                trgt_version = Some(value.clone());
                // Don't break here, continue to find the last occurrence
            }
        }
    }

    // If trgtVersion is not in the header it's either a <v0.5 TRGT VCF or not a TRGT VCF
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
    pub readers: Vec<VcfReader>,
}

impl VcfReaders {
    pub fn new(vcf_files: Vec<PathBuf>) -> Result<Self> {
        let readers = vcf_files
            .into_iter()
            .enumerate()
            .map(|(index, file)| VcfReader::new(file, index))
            .collect::<Result<Vec<_>>>()?;

        Ok(VcfReaders { readers })
    }

    pub fn get_contig_order(&self) -> Result<Vec<String>> {
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

    pub fn merge_headers(&self, dst_header: &mut Header, force_samples: bool) -> Result<()> {
        let mut observed_sample_ids = HashSet::new();

        for reader in &self.readers {
            let src_header = &reader.header;
            unsafe {
                dst_header.inner =
                    rust_htslib::htslib::bcf_hdr_merge(dst_header.inner, src_header.inner);
            }

            for sample_id in src_header.samples() {
                if observed_sample_ids.contains(sample_id) {
                    if force_samples {
                        continue; // If forcing samples, skip duplicates
                    } else {
                        return Err(format!(
                            "Duplicate sample ID found: {}",
                            String::from_utf8_lossy(sample_id)
                        ));
                    }
                }
                observed_sample_ids.insert(sample_id.to_vec());
                dst_header.push_sample(sample_id);
            }
        }

        unsafe {
            rust_htslib::htslib::bcf_hdr_sync(dst_header.inner);
        }

        Ok(())
    }
}
