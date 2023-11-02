use crate::genotype::Ploidy;
use crate::utils::{self, GenomicRegion};
use rust_htslib::faidx;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader, Read as ioRead};
use std::str::FromStr;

#[derive(Debug, PartialEq, Clone)]
pub struct Karyotype {
    ploidy: PloidyInfo,
}

#[derive(Debug, PartialEq, Clone)]
enum PloidyInfo {
    PresetXX,
    PresetXY,
    Custom(HashMap<String, Ploidy>),
}

impl Karyotype {
    pub fn new(encoding: &str) -> Result<Self, String> {
        match encoding {
            "XX" => Ok(Self {
                ploidy: PloidyInfo::PresetXX,
            }),
            "XY" => Ok(Self {
                ploidy: PloidyInfo::PresetXY,
            }),
            _ => Self::from_file(encoding),
        }
    }

    fn from_file(path: &str) -> Result<Self, String> {
        let contents = fs::read_to_string(path).map_err(|e| format!("File {}: {}", path, e))?;

        let ploidies = contents
            .lines()
            .map(|line| {
                let mut parts = line.split_whitespace();
                let chrom = parts.next().ok_or("Missing chromosome".to_string())?;
                let ploidy_str = parts.next().ok_or("Missing ploidy".to_string())?;
                let ploidy = ploidy_str
                    .parse()
                    .map_err(|e: String| format!("Invalid ploidy: {}", e))?;
                Ok((chrom.to_string(), ploidy))
            })
            .collect::<Result<HashMap<_, _>, String>>()?;

        Ok(Self {
            ploidy: PloidyInfo::Custom(ploidies),
        })
    }

    pub fn get_ploidy(&self, chrom: &str) -> Result<Ploidy, String> {
        let is_on_chrx = chrom == "X" || chrom == "chrX";
        let is_on_chry = chrom == "Y" || chrom == "chrY";
        match &self.ploidy {
            PloidyInfo::PresetXX => {
                if is_on_chry {
                    Ok(Ploidy::Zero)
                } else {
                    Ok(Ploidy::Two)
                }
            }
            PloidyInfo::PresetXY => {
                if is_on_chrx || is_on_chry {
                    Ok(Ploidy::One)
                } else {
                    Ok(Ploidy::Two)
                }
            }
            PloidyInfo::Custom(ploidies) => ploidies
                .get(chrom)
                .copied()
                .ok_or_else(|| format!("Ploidy was not specified for chromosome: {}", chrom)),
        }
    }
}

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
    pub region: utils::GenomicRegion,
    pub motifs: Vec<String>,
    pub struc: String,
    pub ploidy: Ploidy,
    pub genotyper: Genotyper,
}

impl Locus {
    pub fn new(
        genome_reader: &faidx::Reader,
        chrom_lookup: &HashSet<String>,
        line: &str,
        flank_len: usize,
        karyotype: &Karyotype,
        genotyper: Genotyper,
    ) -> Result<Self, String> {
        let split_line: Vec<&str> = line.split_whitespace().collect();
        if split_line.len() != 4 {
            return Err(format!("Expected 4 fields, found {}", split_line.len()));
        }

        let (chrom, start, end) = (split_line[0], split_line[1], split_line[2]);
        let region = GenomicRegion::new(&format!("{}:{}-{}", chrom, start, end))?;

        let ploidy = karyotype.get_ploidy(chrom)?;

        let info_fields = split_line[3];
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

        let (left_flank, tr, right_flank) =
            get_tr_and_flanks(genome_reader, chrom_lookup, &region, flank_len)?;

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

pub fn get_loci<'a>(
    catalog_reader: BufReader<Box<dyn ioRead>>,
    genome_reader: &'a faidx::Reader,
    flank_len: usize,
    karyotype: &'a Karyotype,
    genotyper: Genotyper,
) -> impl Iterator<Item = Result<Locus, String>> + 'a {
    let chrom_lookup: HashSet<String> = (0..genome_reader.n_seqs())
        .filter_map(|i| genome_reader.seq_name(i as i32).ok())
        .collect();

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
                    karyotype,
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
    chrom_lookup: &HashSet<String>,
    region: &utils::GenomicRegion,
    flank_len: usize,
) -> Result<(String, String, String), String> {
    let (lf_start, lf_end) = (region.start as usize - flank_len, region.start as usize);
    let (rf_start, rf_end) = (region.end as usize, region.end as usize + flank_len);

    // TODO: This is necessary because faidx is unsafe and segfaults when
    // the region is invalid, so we need to fail gracefully in the event of
    // a bad input. Should be removed if rust_htslib addresses this.
    if !chrom_lookup.contains(&region.contig) {
        return Err(format!(
            "FASTA reference does not contain chromosome '{}' in BED file",
            region.contig
        ));
    }

    let left_flank = genome
        .fetch_seq_string(&region.contig, lf_start, lf_end - 1)
        .map_err(|e| e.to_string())?;
    let tr = genome
        .fetch_seq_string(
            &region.contig,
            region.start as usize,
            region.end as usize - 1,
        )
        .map_err(|e| e.to_string())?;
    let right_flank = genome
        .fetch_seq_string(&region.contig, rf_start, rf_end - 1)
        .map_err(|e| e.to_string())?;
    Ok((
        left_flank.to_uppercase(),
        tr.to_uppercase(),
        right_flank.to_uppercase(),
    ))
}

fn decode_fields(info_fields: &str) -> Result<HashMap<&str, String>, String> {
    let mut fields = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding)?;
        if fields.insert(name, value.to_string()).is_some() {
            return Err(format!("Duplicate field: {}", name));
        }
    }
    Ok(fields)
}

fn decode_info_field(encoding: &str) -> Result<(&str, &str), String> {
    if encoding.is_empty() {
        return Err("Field is empty".to_string());
    }

    let mut name_and_value = encoding.splitn(2, '=');
    let error_message = || format!("Invalid field entry: {}", encoding);

    let name = name_and_value
        .next()
        .filter(|s| !s.is_empty())
        .ok_or_else(error_message)?;

    let value = name_and_value
        .next()
        .filter(|s| !s.is_empty())
        .ok_or_else(error_message)?;

    Ok((name, value))
}
