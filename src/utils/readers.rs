use super::Result;
use flate2::read::MultiGzDecoder;
use rust_htslib::faidx;
use std::fs::File;
use std::io::{BufReader, Read as ioRead};
use std::path::Path;

pub fn open_catalog_reader(path: &Path) -> Result<BufReader<Box<dyn ioRead>>> {
    fn is_gzipped(path: &Path) -> bool {
        let path_str = path.to_string_lossy().to_lowercase();
        path_str.ends_with(".gz") || path_str.ends_with(".gzip")
    }
    let file = File::open(path).map_err(|e| e.to_string())?;
    if is_gzipped(path) {
        let gz_decoder = MultiGzDecoder::new(file);
        if gz_decoder.header().is_some() {
            Ok(BufReader::new(Box::new(gz_decoder)))
        } else {
            Err(format!("Invalid gzip header: {}", path.to_string_lossy()))
        }
    } else {
        Ok(BufReader::new(Box::new(file)))
    }
}

pub fn open_genome_reader(path: &Path) -> Result<faidx::Reader> {
    let extension = path.extension().unwrap().to_str().unwrap();
    let fai_path = path.with_extension(extension.to_owned() + ".fai");
    if !fai_path.exists() {
        return Err(format!(
            "Reference index file not found: {}. Create it using 'samtools faidx {}'",
            fai_path.display(),
            path.display()
        ));
    }
    faidx::Reader::from_path(path).map_err(|e| e.to_string())
}
