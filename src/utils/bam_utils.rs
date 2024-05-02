use crate::utils::Result;
use rust_htslib::bam::{self, Read};
use std::{collections::HashSet, path::Path};

pub fn get_bam_header(bam_path: &Path) -> Result<bam::Header> {
    let bam = bam::IndexedReader::from_path(bam_path)
        .map_err(|e| format!("Failed to create bam reader: {}", e))?;
    Ok(bam::Header::from_template(bam.header()))
}

pub fn is_bam_mapped(bam_header: &bam::Header) -> bool {
    // input is already sorted because it fails an index.
    // If it is mapped, the index needs the SQ tags to fetch data.
    for line in String::from_utf8(bam_header.to_bytes()).unwrap().lines() {
        if line.starts_with("@SQ") {
            return true;
        }
    }
    false
}

pub fn get_sample_name(reads_path: &Path, bam_header: &bam::Header) -> Result<String> {
    let header_hashmap = bam_header.to_hashmap();
    let mut sample_names = HashSet::new();

    if let Some(rg_fields) = header_hashmap.get("RG") {
        for rg_field in rg_fields {
            if let Some(sample_name) = rg_field.get("SM") {
                sample_names.insert(sample_name.to_owned());
            }
        }
    }

    match sample_names.len() {
        1 => return Ok(sample_names.into_iter().next().unwrap()),
        0 => log::warn!("No sample names found"),
        _ => log::warn!("Multiple sample names found"),
    };

    let sample = reads_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("Invalid reads file name")?
        .to_string();

    Ok(sample)
}
