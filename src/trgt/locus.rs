use crate::trgt::{locus_group::LocusGroup, reads::LocusRead};
use crate::utils::{
    open_catalog_reader, open_genome_reader, GenomicRegion, Genotyper, InputSource, Karyotype,
    Ploidy, Result,
};
use crossbeam_channel::Sender;
use rust_htslib::faidx;
use std::{collections::HashMap, io::BufRead};

#[derive(Debug)]
pub struct Locus {
    pub id: String,
    pub left_flank: Vec<u8>,
    pub tr: Vec<u8>,
    pub right_flank: Vec<u8>,
    pub region: GenomicRegion,
    pub motifs: Vec<Vec<u8>>,
    pub struc: String,
    pub ploidy: Ploidy,
    pub genotyper: Genotyper,
    pub reads: Vec<LocusRead>,
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

        let (chrom, _start, _end, info_fields) = match &split_line[..] {
            [chrom, start, end, info_fields] => (*chrom, *start, *end, *info_fields),
            _ => unreachable!(),
        };

        let region = GenomicRegion::from_bed_fields(&split_line[..3])?;

        check_region_bounds(&region, flank_len, chrom_lookup)?;

        let ploidy = karyotype.get_ploidy(chrom)?;
        let fields = decode_fields(info_fields)?;
        let id = get_field(&fields, "ID")?;
        let motifs = get_field(&fields, "MOTIFS")?
            .split(',')
            .map(|s| s.as_bytes().to_vec())
            .collect();
        // TODO: This should not be mandatory anymore
        let struc = get_field(&fields, "STRUC")?;

        let (left_flank, tr, right_flank) = get_tr_and_flanks(genome_reader, &region, flank_len)?;

        Ok(Locus {
            id,
            left_flank,
            tr,
            right_flank,
            region,
            motifs,
            ploidy,
            struc,
            genotyper,
            reads: Vec::new(),
        })
    }

    #[cfg(test)]
    pub fn test_base() -> Self {
        Self {
            id: "locus".into(),
            left_flank: b"A".to_vec(),
            tr: b"T".to_vec(),
            right_flank: b"C".to_vec(),
            region: GenomicRegion {
                contig: "chr1".into(),
                start: 1,
                end: 2,
            },
            motifs: vec![b"T".to_vec()],
            struc: "1".into(),
            ploidy: Ploidy::Two,
            genotyper: Genotyper::Size,
            reads: Vec::new(),
        }
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

#[allow(clippy::too_many_arguments)]
pub fn stream_locus_groups_into_channel(
    repeats_src: &InputSource,
    genome_src: &InputSource,
    flank_len: usize,
    genotyper: Genotyper,
    karyotype: &Karyotype,
    locus_group_tx: Sender<Result<LocusGroup>>,
    max_group_span: u32,
    max_group_size: usize,
) -> Result<()> {
    let mut catalog_reader = open_catalog_reader(repeats_src).map_err(|e| e.to_string())?;
    let genome_reader = open_genome_reader(genome_src).map_err(|e| e.to_string())?;
    let chrom_lookup = create_chrom_lookup(&genome_reader).map_err(|e| e.to_string())?;

    let mut current_group: Option<LocusGroup> = None;
    let mut line = String::new();
    let mut line_number = 0;

    while {
        line.clear();
        match catalog_reader.read_line(&mut line) {
            Ok(0) => false, // EOF
            Ok(_) => true,
            Err(e) => {
                if !send_locus_group_result(
                    &locus_group_tx,
                    Err(format!("Error reading BED line {}: {}", line_number + 1, e)),
                ) {
                    return Ok(());
                }
                false
            }
        }
    } {
        line_number += 1;

        let locus_result = Locus::new(
            &genome_reader,
            &chrom_lookup,
            &line,
            flank_len,
            karyotype,
            genotyper,
        );

        let locus = match locus_result {
            Ok(locus) => locus,
            Err(e) => {
                if !send_locus_group_result(
                    &locus_group_tx,
                    Err(format!("Error at BED line {}: {}", line_number, e)),
                ) {
                    return Ok(());
                }
                continue;
            }
        };

        if let Some(group) = &mut current_group {
            if group.len() < max_group_size && group.can_add_locus(&locus, max_group_span) {
                group.add_locus_unchecked(locus);
                continue;
            }
            if !send_locus_group_result(&locus_group_tx, Ok(current_group.take().unwrap())) {
                return Ok(());
            }
        }
        current_group = Some(LocusGroup::new(locus));
    }

    if let Some(group) = current_group.take() {
        send_locus_group_result(&locus_group_tx, Ok(group));
    }

    Ok(())
}

#[inline]
fn send_locus_group_result(
    sender: &Sender<Result<LocusGroup>>,
    result: Result<LocusGroup>,
) -> bool {
    match sender.send(result) {
        Ok(()) => true,
        Err(_) => {
            log::error!("Stopped locus grouping because downstream receiver disconnected.");
            false
        }
    }
}

pub fn get_tr_and_flanks(
    genome: &faidx::Reader,
    region: &GenomicRegion,
    flank_len: usize,
) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let full_start = region.start as usize - flank_len;
    let full_end = region.end as usize + flank_len - 1;

    let mut full_seq = genome
        .fetch_seq(&region.contig, full_start, full_end)
        .map_err(|e| {
            format!(
                "Error fetching sequence for region {}:{}-{}: {}",
                &region.contig, full_start, full_end, e
            )
        })?;

    full_seq.make_ascii_uppercase();

    let tr_start = flank_len;
    let tr_end = flank_len + (region.end - region.start) as usize;
    let right_flank_start = tr_end;

    let left_flank = full_seq[..tr_start].to_vec();
    let tr = full_seq[tr_start..tr_end].to_vec();
    let right_flank = full_seq[right_flank_start..].to_vec();

    Ok((left_flank, tr, right_flank))
}

pub fn get_field(fields: &HashMap<&str, &str>, key: &str) -> Result<String> {
    fields
        .get(key)
        .ok_or_else(|| format!("{} field missing", key))
        .map(|s| s.to_string())
}

pub fn decode_fields(info_fields: &str) -> Result<HashMap<&str, &str>> {
    let mut fields = HashMap::new();
    for field_encoding in info_fields.split(';') {
        let (name, value) = decode_info_field(field_encoding)?;
        if fields.insert(name, value).is_some() {
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

pub fn check_region_bounds(
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_fields_accepts_unrecognized_catalog_info_fields() {
        let fields = decode_fields("ID=TR1;GENE=HTT;MOTIFS=CAG;SOURCE=custom;STRUC=<TR>").unwrap();
        assert_eq!(get_field(&fields, "ID").unwrap(), "TR1");
        assert_eq!(get_field(&fields, "MOTIFS").unwrap(), "CAG");
        assert_eq!(get_field(&fields, "STRUC").unwrap(), "<TR>");
        assert_eq!(fields.get("GENE"), Some(&"HTT"));
        assert_eq!(fields.get("SOURCE"), Some(&"custom"));
    }
}
