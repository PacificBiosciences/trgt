use crate::utils::Result;
use rust_htslib::{
    bam::{self, Read},
    bcf,
};
use std::path::Path;

pub fn create_writer<T, F>(output_prefix: &Path, output_suffix: &str, f: F) -> Result<T>
where
    F: FnOnce(&str) -> Result<T>,
{
    let mut output_path = output_prefix.to_path_buf().into_os_string();
    output_path.push(format!(".{output_suffix}"));
    f(output_path.to_str().unwrap())
}

pub fn open_vcf_reader(path: &Path) -> Result<bcf::Reader> {
    let vcf = match bcf::Reader::from_path(path) {
        Ok(vcf) => vcf,
        Err(e) => return Err(format!("Failed to open VCF file {}: {}", path.display(), e)),
    };
    Ok(vcf)
}

pub fn open_bam_reader(path: &Path, threads: usize) -> Result<bam::IndexedReader> {
    let mut bam = match bam::IndexedReader::from_path(path) {
        Ok(bam) => bam,
        Err(e) => return Err(format!("Failed to open BAM file {}: {}", path.display(), e)),
    };
    bam.set_threads(threads).map_err(|e| e.to_string())?;
    Ok(bam)
}
