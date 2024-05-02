use crate::utils::{
    open_catalog_reader, open_genome_reader, GenomicRegion, Genotyper, Karyotype, Ploidy, Result,
};
use crossbeam_channel::Sender;
use rust_htslib::faidx;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader, Read as ioRead},
    path::Path,
};

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

pub fn create_chrom_lookup(reader: &faidx::Reader) -> Result<HashMap<String, u32>> {
    let num_seqs = reader.n_seqs() as usize;
    let mut map = HashMap::with_capacity(num_seqs);
    for i in 0..num_seqs {
        let name = reader.seq_name(i as i32).map_err(|e| e.to_string())?;
        let len = reader.fetch_seq_len(&name);
        let len_u32 = u32::try_from(len).map_err(|_| {
            format!(
                "Sequence length for '{}' is negative and cannot be converted to u32",
                &name
            )
        })?;
        map.insert(name, len_u32);
    }
    Ok(map)
}

pub fn stream_loci_into_channel(
    repeats_path: &Path,
    genome_path: &Path,
    flank_len: usize,
    genotyper: Genotyper,
    karyotype: &Karyotype,
    sender: Sender<Result<Locus>>,
) {
    let catalog_reader = open_catalog_reader(repeats_path).unwrap();
    let genome_reader = open_genome_reader(genome_path).unwrap();
    let chrom_lookup = create_chrom_lookup(&genome_reader).unwrap();

    for (line_number, result_line) in catalog_reader.lines().enumerate() {
        let line = match result_line {
            Ok(line) => line,
            Err(err) => {
                let error = format!("Error at BED line {}: {}", line_number + 1, err);
                sender
                    .send(Err(error))
                    .expect("Failed to send error through channel");
                return;
            }
        };

        let locus_result = Locus::new(
            &genome_reader,
            &chrom_lookup,
            &line,
            flank_len,
            karyotype,
            genotyper,
        );
        let locus = match locus_result {
            Ok(locus) => Ok(locus),
            Err(e) => {
                let error = format!("Error at BED line {}: {}", line_number + 1, e);
                sender
                    .send(Err(error))
                    .expect("Failed to send error through channel");
                continue;
            }
        };

        sender
            .send(locus)
            .expect("Failed to send locus through channel");
    }
}

pub fn get_loci(
    catalog_reader: BufReader<Box<dyn ioRead>>,
    genome_reader: faidx::Reader,
    karyotype: Karyotype,
    flank_len: usize,
    genotyper: Genotyper,
) -> impl Iterator<Item = Result<Locus>> {
    let chrom_lookup = create_chrom_lookup(&genome_reader).unwrap();

    catalog_reader
        .lines()
        .enumerate()
        .map(move |(line_number, result_line)| {
            result_line
                .map_err(|e| format!("Error at BED line {}: {}", line_number + 1, e))
                .and_then(|line| {
                    Locus::new(
                        &genome_reader,
                        &chrom_lookup,
                        &line,
                        flank_len,
                        &karyotype,
                        genotyper,
                    )
                    .map_err(|e| format!("Error at BED line {}: {}", line_number + 1, e))
                })
        })
}

fn get_tr_and_flanks(
    genome: &faidx::Reader,
    region: &GenomicRegion,
    flank_len: usize,
) -> Result<(String, String, String)> {
    let fetch_seq = |start: usize, end: usize| {
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

    let left_flank = fetch_seq(region.start as usize - flank_len, region.start as usize - 1)?;
    let tr = fetch_seq(region.start as usize, region.end as usize - 1)?;
    let right_flank = fetch_seq(region.end as usize, region.end as usize + flank_len - 1)?;

    Ok((left_flank, tr, right_flank))
}

fn decode_fields(info_fields: &str) -> Result<HashMap<&str, String>> {
    let mut fields = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding).map_err(|e| e.to_string())?;
        if fields.insert(name, value.to_string()).is_some() {
            return Err(format!("Duplicate field name: '{}'", name));
        }
    }
    Ok(fields)
}

fn decode_info_field(encoding: &str) -> Result<(&str, &str)> {
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
) -> Result<()> {
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
