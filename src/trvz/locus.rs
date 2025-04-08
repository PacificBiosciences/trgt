use crate::{
    trgt::locus::{check_region_bounds, decode_fields, get_field, get_tr_and_flanks},
    trvz::struc::RegionLabel,
    utils::{GenomicRegion, Result},
};
use rust_htslib::faidx;
use std::collections::HashMap;

#[derive(Debug)]
pub struct Allele {
    pub seq: String,
    pub region_labels: Vec<RegionLabel>,
    pub flank_labels: Vec<RegionLabel>,
    //pub base_labels: Vec<BaseLabel>,
}

#[derive(Debug)]
pub struct Locus {
    pub id: String,
    pub struc: String,
    pub motifs: Vec<String>,
    pub left_flank: String,
    pub right_flank: String,
    pub region: GenomicRegion,
}

impl Locus {
    pub fn new(
        genome_reader: &faidx::Reader,
        chrom_lookup: &HashMap<String, u32>,
        line: &str,
        flank_len: usize,
    ) -> Result<Self> {
        const EXPECTED_FIELD_COUNT: usize = 4;
        let split_line: Vec<&str> = line.split_whitespace().collect();
        if split_line.len() != EXPECTED_FIELD_COUNT {
            return Err(format!(
                "Expected {} fields in the format 'chrom start end info', found {}: {}",
                EXPECTED_FIELD_COUNT,
                split_line.len(),
                line
            ));
        }

        let (chrom, start, end, info_fields) = match &split_line[..] {
            [chrom, start, end, info_fields] => (*chrom, *start, *end, *info_fields),
            _ => unreachable!(),
        };

        let region = GenomicRegion::from_string(&format!("{}:{}-{}", chrom, start, end))?;
        check_region_bounds(&region, flank_len, chrom_lookup)?;

        let fields = decode_fields(info_fields)?;
        let id = get_field(&fields, "ID")?;
        let motifs = get_field(&fields, "MOTIFS")?
            .split(',')
            .map(|s| s.to_string())
            .collect();
        let struc = get_field(&fields, "STRUC")?;

        let (left_flank, _, right_flank) = get_tr_and_flanks(genome_reader, &region, flank_len)?;

        Ok(Locus {
            id,
            struc,
            motifs,
            left_flank,
            right_flank,
            region,
        })
    }
}
