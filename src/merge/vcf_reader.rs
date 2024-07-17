use crate::utils::Result;
use rust_htslib::bcf::{self, header::HeaderView, Header, HeaderRecord, Read};
use semver::Version;
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

pub struct VcfReader {
    pub reader: bcf::IndexedReader,
    pub header: bcf::header::HeaderView,
    pub current_record: bcf::Record,
    pub version: Version,
    pub index: usize,
}

impl VcfReader {
    pub fn new(file: PathBuf, index: usize) -> Result<Self> {
        log::info!("Start loading VCF {:?}", &file);
        // TODO: Check if file is a VCF
        // TODO: Check if indexed VCF
        // TODO: Check if valid VCF
        let reader = bcf::IndexedReader::from_path(&file)
            .map_err(|e| format!("Failed to open VCF file {}: {}", file.display(), e))?;
        let header = reader.header().clone();

        let version = get_trgt_version(&header, &file)?;
        log::debug!("{:?} has version: {}", file.file_name().unwrap(), version);

        if header.sample_count() > 1 {
            return Err(format!(
                "Unsupported: VCF file with multiple samples: {}",
                file.display()
            ));
        }

        // TODO: Create a normalized struct for variant records
        let current_record = reader.empty_record();
        log::info!("Finished loading VCF {:?}", &file);
        Ok(VcfReader {
            reader,
            header,
            current_record,
            version,
            index,
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
    // TODO: Add logic to deal with merged TRGT VCFs (assume latest version?)
    let mut trgt_version = None;

    for record in vcf_header.header_records().iter() {
        if let HeaderRecord::Generic { key, value } = record {
            if key == "trgtVersion" {
                trgt_version = Some(value.clone());
                break;
            }
        }
    }

    // If trgtVersion is not in the header its either a <v0.5 TRGT VCF or not a TRGT VCF
    if trgt_version.is_none() {
        let mut has_allr = false;
        let mut has_alci = false;
        let mut is_integer_am = false;
        for record in vcf_header.header_records().iter() {
            if let HeaderRecord::Format { key: _, values } = record {
                if let Some(id) = values.get("ID") {
                    if id == "ALLR" {
                        has_allr = true;
                    }
                    if id == "ALCI" {
                        has_alci = true;
                    }
                    if id == "AM" {
                        is_integer_am = values.get("Type").unwrap() == "Integer";
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

    let version = Version::parse(&trgt_version.unwrap()).unwrap();

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

    pub fn merge_headers(&self, dst_header: &mut Header, force_samples: bool) -> Result<()> {
        let mut observed_sample_ids = HashSet::new();

        for reader in &self.readers {
            let src_header = &reader.header;
            // TODO: error handling
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
