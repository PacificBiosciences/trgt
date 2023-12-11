use crate::faidx;
use crate::genotype::Ploidy;
use crate::karyotype::Karyotype;
use crate::utils::GenomicRegion;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read as ioRead};
use std::str::FromStr;

#[derive(Debug, Clone, Copy)]
pub enum Genotyper {
    Size,
    Cluster,
}

impl FromStr for Genotyper {
    type Err = &'static str;
    fn from_str(genotyper: &str) -> Result<Self, Self::Err> {
        match genotyper {
            "size" => Ok(Genotyper::Size),
            "cluster" => Ok(Genotyper::Cluster),
            _ => Err("Invalid genotyper"),
        }
    }
}

#[derive(Debug)]
pub struct Locus {
    pub id: String,
    pub left_flank: String,
    pub tr: String,
    pub right_flank: String,
    pub region: GenomicRegion,
    pub motifs: Vec<String>,
    pub struc: String,
    pub ploidy: Ploidy,
    pub genotyper: Genotyper,
}

impl Locus {
    pub fn new(
        genome_reader: &faidx::Reader,
        chrom_lookup: &HashMap<String, u32>,
        line: &str,
        flank_len: usize,
        karyotype: &Karyotype,
        genotyper: Genotyper,
    ) -> Result<Self, String> {
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

        let region = GenomicRegion::new(&format!("{}:{}-{}", chrom, start, end))?;

        check_region_bounds(&region, flank_len, chrom_lookup)?;

        let ploidy = karyotype.get_ploidy(chrom)?;

        let fields = decode_fields(info_fields)?;

        let get_field = |key: &str| {
            fields
                .get(key)
                .ok_or_else(|| format!("{} field missing", key))
                .map(|s| s.to_string())
        };

        let id = get_field("ID")?;
        let motifs = get_field("MOTIFS")?
            .split(',')
            .map(|s| s.to_string())
            .collect();
        let struc = get_field("STRUC")?;

        let (left_flank, tr, right_flank) = get_tr_and_flanks(genome_reader, &region, flank_len)?;

        Ok(Locus {
            id,
            left_flank,
            tr,
            right_flank,
            region,
            motifs,
            struc,
            ploidy,
            genotyper,
        })
    }
}

pub fn get_loci(
    catalog_reader: BufReader<Box<dyn ioRead>>,
    genome_reader: &faidx::Reader,
    karyotype: Karyotype,
    flank_len: usize,
    genotyper: Genotyper,
) -> impl Iterator<Item = Result<Locus, String>> + '_ {
    let chrom_lookup = genome_reader.create_chrom_lookup().unwrap();

    catalog_reader
        .lines()
        .enumerate()
        .filter_map(move |(line_number, result_line)| match result_line {
            Ok(line) => {
                match Locus::new(
                    genome_reader,
                    &chrom_lookup,
                    &line,
                    flank_len,
                    &karyotype,
                    genotyper,
                ) {
                    Ok(locus) => Some(Ok(locus)),
                    Err(e) => Some(Err(format!("Error at BED line {}: {}", line_number + 1, e))),
                }
            }
            Err(e) => Some(Err(format!("Error at BED line {}: {}", line_number + 1, e))),
        })
}

fn get_tr_and_flanks(
    genome: &faidx::Reader,
    region: &GenomicRegion,
    flank_len: usize,
) -> Result<(String, String, String), String> {
    let fetch_flank = |start: usize, end: usize| {
        genome
            .fetch_seq_string(&region.contig, start, end)
            .map_err(|e| {
                format!(
                    "Error fetching sequence for region {}:{}-{}: {}",
                    &region.contig, start, end, e
                )
            })
            .map(|seq| seq.to_uppercase())
    };

    let left_flank = fetch_flank(region.start as usize - flank_len, region.start as usize - 1)?;
    let tr = fetch_flank(region.start as usize, region.end as usize - 1)?;
    let right_flank = fetch_flank(region.end as usize, region.end as usize + flank_len - 1)?;

    Ok((left_flank, tr, right_flank))
}

fn decode_fields(info_fields: &str) -> Result<HashMap<&str, String>, String> {
    let mut fields = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding).map_err(|e| e.to_string())?;
        if fields.insert(name, value.to_string()).is_some() {
            return Err(format!("Duplicate field name: '{}'", name));
        }
    }
    Ok(fields)
}

fn decode_info_field(encoding: &str) -> Result<(&str, &str), String> {
    let error_message = || format!("Field must be in 'name=value' format: '{}'", encoding);
    let parts: Vec<&str> = encoding.splitn(2, '=').collect();
    if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
        Err(error_message())
    } else {
        Ok((parts[0], parts[1]))
    }
}

fn check_region_bounds(
    region: &GenomicRegion,
    flank_len: usize,
    chrom_lookup: &HashMap<String, u32>,
) -> Result<(), String> {
    let chrom_length = *chrom_lookup.get(&region.contig).ok_or_else(|| {
        format!(
            "FASTA reference does not contain chromosome '{}' in BED file",
            &region.contig
        )
    })?;

    let flank_len_u32 = flank_len as u32;

    if region.start < flank_len_u32 + 1 {
        return Err(format!(
            "Region start '{}' with flank length '{}' underflows for chromosome '{}'.",
            region.start, flank_len, &region.contig
        ));
    }

    let adjusted_end = region.end.checked_add(flank_len_u32).ok_or_else(|| {
        format!(
            "Region end '{}' with flank length '{}' overflows for chromosome '{}'.",
            region.end, flank_len, &region.contig
        )
    })?;

    if adjusted_end > chrom_length {
        return Err(format!(
            "Region end '{}' with flank length '{}' exceeds chromosome '{}' bounds (0..{}).",
            adjusted_end, flank_len, &region.contig, chrom_length
        ));
    }

    Ok(())
}
