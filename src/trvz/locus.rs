use crate::trvz::struc::RegionLabel;
use crate::utils::{GenomicRegion, Result};
use itertools::Itertools;
use rust_htslib::faidx;
use std::collections::HashMap;

#[derive(Debug)]
pub struct Locus {
    pub id: String,
    pub struc: String,
    pub motifs: Vec<String>,
    pub left_flank: String,
    pub right_flank: String,
    pub region: GenomicRegion,
}

#[derive(Debug)]
pub struct Allele {
    pub seq: String,
    pub region_labels: Vec<RegionLabel>,
    pub flank_labels: Vec<RegionLabel>,
    //pub base_labels: Vec<BaseLabel>,
}

fn get_flanks(
    genome: &faidx::Reader,
    region: &GenomicRegion,
    flank_len: usize,
) -> Result<(String, String)> {
    let (lf_start, lf_end) = (region.start as usize - flank_len, region.start as usize);
    let (rf_start, rf_end) = (region.end as usize, region.end as usize + flank_len);

    let left_flank = match genome.fetch_seq_string(&region.contig, lf_start, lf_end - 1) {
        Ok(seq) => seq,
        Err(_) => return Err(format!("Unable to extract: {:?}", region)),
    };
    let right_flank = match genome.fetch_seq_string(&region.contig, rf_start, rf_end - 1) {
        Ok(seq) => seq,
        Err(_) => return Err(format!("Unable to extract: {:?}", region)),
    };
    Ok((left_flank.to_uppercase(), right_flank.to_uppercase()))
}

fn decode_info_field(encoding: &str) -> Result<(&str, &str)> {
    let name_and_value: Vec<&str> = encoding.split('=').collect();
    if name_and_value.len() != 2 {
        return Err(format!("Invalid entry: {encoding}"));
    }
    Ok((name_and_value[0], name_and_value[1]))
}

pub fn decode(flank_len: usize, genome: &faidx::Reader, encoding: &str) -> Result<Locus> {
    let split_line: Vec<&str> = encoding.split_whitespace().collect();
    if split_line.len() != 4 {
        return Err(format!("Invalid bed entry encountered: {encoding}"));
    }

    let encoding = format!("{}:{}-{}", split_line[0], split_line[1], split_line[2]);
    let region = GenomicRegion::from_string(&encoding).unwrap();

    let info_fields = split_line[3];

    let mut name_to_value = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding)?;
        name_to_value.insert(
            name,
            match name {
                "ID" => value,
                "STRUC" => value,
                "MOTIFS" => value,
                _ => "",
            },
        );
    }

    let id = name_to_value
        .get("ID")
        .ok_or_else(|| format!("Missing ID in {}", encoding))?
        .to_string();

    let motifs = name_to_value
        .get("MOTIFS")
        .ok_or_else(|| format!("Missing MOTIFS in {}", encoding))?
        .split(',')
        .map(|s| s.to_string())
        .collect_vec();

    let struc = name_to_value
        .get("STRUC")
        .ok_or_else(|| format!("Missing STRUC in {}", encoding))?
        .to_string();

    let (left_flank, right_flank) = get_flanks(genome, &region, flank_len)?;

    Ok(Locus {
        id,
        struc,
        motifs,
        left_flank,
        right_flank,
        region,
    })
}
