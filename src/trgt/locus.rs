use crate::trgt::{locus_group::LocusGroup, reads::LocusRead};
use crate::utils::{
    fnv1a64, open_catalog_reader, open_genome_reader, GenomicRegion, Genotyper, InputSource,
    Karyotype, Ploidy, Result,
};
use crossbeam_channel::Sender;
use rust_htslib::faidx;
use std::{collections::HashMap, io::BufRead};

/// How many records a previous run left for one locus ID, and how many of them are still waiting to
/// be matched to a line of the catalog.
///
/// `digest` combines the definition of every record kept for the ID with the definition of every
/// catalog line matched to one of them. The two sides cancel each other out, so a digest of zero at
/// the end of the catalog means the catalog still describes those loci the way the records do.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct KeptRecords {
    pub total: u32,
    pub unmatched: u32,
    pub digest: u64,
}

impl KeptRecords {
    /// Counts one more record kept for the locus ID, defined as the given locus.
    pub fn add_record(&mut self, definition: &LocusDefinition) {
        self.total += 1;
        self.unmatched = self.total;
        self.digest ^= definition.digest();
    }
}

/// What a catalog says a locus is: where it lies and what it is built of.
///
/// Genotypes follow all of it, so a resumed run keeps the records of a previous run only while the
/// catalog still defines their loci this way.
#[derive(Debug, Clone, PartialEq)]
pub struct LocusDefinition {
    pub region: GenomicRegion,
    pub motifs: String,
    pub struc: String,
}

impl LocusDefinition {
    /// Returns a value standing for the definition, to recognize it by later.
    ///
    /// Digests are combined with exclusive or, which ignores the order the loci come in. That is
    /// what lets the records of a previous run, written in the order genotyping finished them, be
    /// compared against the catalog, which is in its own order.
    fn digest(&self) -> u64 {
        fnv1a64(
            format!(
                "{}:{}-{}:{}:{}",
                self.region.contig, self.region.start, self.region.end, self.motifs, self.struc
            )
            .as_bytes(),
        )
    }
}

/// Returns what a catalog BED line defines its locus as, or `None` if the line is not a locus.
fn parse_locus_definition(line: &str) -> Option<LocusDefinition> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    Some(LocusDefinition {
        region: GenomicRegion::from_bed_fields(fields.get(..3)?).ok()?,
        motifs: parse_info_field(line, "MOTIFS")?.to_string(),
        struc: parse_info_field(line, "STRUC")?.to_string(),
    })
}

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
    mut completed_ids: HashMap<String, KeptRecords>,
) -> Result<()> {
    let mut catalog_reader = open_catalog_reader(repeats_src).map_err(|e| e.to_string())?;
    let genome_reader = open_genome_reader(genome_src).map_err(|e| e.to_string())?;
    let chrom_lookup = create_chrom_lookup(&genome_reader).map_err(|e| e.to_string())?;

    let mut current_group: Option<LocusGroup> = None;
    let mut line = String::new();
    let mut line_number = 0;
    let mut num_skipped = 0;

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

        // Skip a locus once for every record the previous run left for its ID.
        if !completed_ids.is_empty() {
            let is_completed = take_completed_record(&mut completed_ids, &line)
                .map_err(|e| format!("Error at BED line {line_number}: {e}"))?;
            if is_completed {
                num_skipped += 1;
                continue;
            }
        }

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

    if num_skipped > 0 {
        log::info!("Skipped {num_skipped} of the {line_number} catalog loci that the previous run already genotyped");
    }

    Ok(())
}

/// Checks that every locus of the catalog that a previous run has a record for can be told apart
/// from the others.
///
/// A record names its locus by ID alone, and results are written in the order genotyping finishes
/// them rather than in catalog order, so when the catalog gives an ID to more loci than the
/// previous run left records for, there is no way to tell which of those loci still need
/// genotyping. Reading the catalog for this before the previous run's output is replaced is what
/// lets such a catalog stop the run with the output file untouched.
pub fn check_catalog_ids(
    repeats_src: &InputSource,
    completed_ids: &mut HashMap<String, KeptRecords>,
) -> Result<()> {
    if completed_ids.is_empty() {
        return Ok(());
    }

    let mut catalog_reader = open_catalog_reader(repeats_src).map_err(|e| e.to_string())?;
    let mut line = String::new();
    let mut line_number = 0;
    loop {
        line.clear();
        match catalog_reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => line_number += 1,
            Err(e) => return Err(format!("Error reading BED line {}: {}", line_number + 1, e)),
        }

        let Some(id) = parse_locus_id(&line) else {
            continue;
        };
        let definition = parse_locus_definition(&line);
        let Some(kept) = completed_ids.get_mut(id) else {
            continue;
        };
        if kept.unmatched == 0 {
            return Err(format!(
                "Error at BED line {line_number}: {}",
                reused_id_error(id)
            ));
        }
        kept.unmatched -= 1;
        if let Some(definition) = definition {
            kept.digest ^= definition.digest();
        }
    }

    for (id, kept) in completed_ids.iter_mut() {
        if kept.unmatched != 0 || kept.digest != 0 {
            return Err(format!(
                "cannot resume a run whose catalog no longer defines locus '{id}' the way the \
                 output of the previous run does, because that locus has been moved, resized, given \
                 other motifs or another structure, or taken out of the catalog since. Resume with \
                 the catalog the run started with, or rerun without --resume"
            ));
        }
        kept.unmatched = kept.total;
    }
    Ok(())
}

/// Accounts for one record that a previous run already wrote for the locus of a catalog BED line.
///
/// Returns whether the line has an ID that the previous run still has an unaccounted record for, in
/// which case that record is counted off.
///
/// Returns an error when an earlier line already used up the records of this line's ID, which
/// [`check_catalog_ids`] only misses if the catalog changed since it read it.
fn take_completed_record(
    completed_ids: &mut HashMap<String, KeptRecords>,
    line: &str,
) -> Result<bool> {
    let Some(id) = parse_locus_id(line) else {
        return Ok(false);
    };
    let Some(kept) = completed_ids.get_mut(id) else {
        return Ok(false);
    };
    if kept.unmatched == 0 {
        return Err(reused_id_error(id));
    }
    kept.unmatched -= 1;
    Ok(true)
}

fn reused_id_error(id: &str) -> String {
    format!(
        "cannot resume a run whose catalog gives the ID '{id}' to more loci than the output of the \
         previous run has records for, because nothing says which of those loci it already \
         genotyped. Repeat IDs must be unique. Give the loci distinct IDs, or rerun without --resume"
    )
}

/// Returns the value of the ID field of a catalog BED line, or `None` if the line has no ID field.
fn parse_locus_id(line: &str) -> Option<&str> {
    parse_info_field(line, "ID")
}

/// Returns the value of one of the info fields of a catalog BED line.
fn parse_info_field<'a>(line: &'a str, name: &str) -> Option<&'a str> {
    line.split_whitespace()
        .nth(3)?
        .split(';')
        .find_map(|field| field.strip_prefix(name)?.strip_prefix('='))
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
    use std::str::FromStr;

    /// The records a previous run left for one ID, none of them matched to a catalog line yet.
    fn kept_records(total: u32) -> KeptRecords {
        KeptRecords {
            total,
            unmatched: total,
            digest: 0,
        }
    }

    /// The record a previous run left for a locus defined the way `write_catalog` defines them.
    fn kept_record_at(contig: &str, start: u32, end: u32) -> KeptRecords {
        let mut kept = KeptRecords::default();
        kept.add_record(&test_definition(contig, start, end));
        kept
    }

    /// A locus defined the way `write_catalog` defines them.
    fn test_definition(contig: &str, start: u32, end: u32) -> LocusDefinition {
        LocusDefinition {
            region: GenomicRegion {
                contig: contig.to_string(),
                start,
                end,
            },
            motifs: "CAG".to_string(),
            struc: "<TR>".to_string(),
        }
    }

    /// The region `write_catalog` gives the locus at `index`.
    fn catalog_region(index: u32) -> (u32, u32) {
        let start = 100 + index * 100;
        (start, start + 30)
    }

    /// Writes a catalog and returns it as an input source, one line per given ID.
    fn write_catalog(dir: &std::path::Path, ids: &[&str]) -> InputSource {
        let path = dir.join("repeats.bed");
        let lines: String = ids
            .iter()
            .enumerate()
            .map(|(index, id)| {
                let start = 100 + index * 100;
                format!(
                    "chr1\t{start}\t{}\tID={id};MOTIFS=CAG;STRUC=<TR>\n",
                    start + 30
                )
            })
            .collect();
        std::fs::write(&path, lines).unwrap();
        InputSource::from_str(path.to_str().unwrap()).unwrap()
    }

    #[test]
    fn take_completed_record_skips_a_locus_once_per_record_of_the_previous_run() {
        let line = "chr1\t100\t200\tID=TR1;MOTIFS=CAG;STRUC=<TR>\n";
        let mut completed_ids = HashMap::from([("TR1".to_string(), kept_records(2))]);

        assert!(take_completed_record(&mut completed_ids, line).unwrap());
        assert!(take_completed_record(&mut completed_ids, line).unwrap());
        assert_eq!(completed_ids["TR1"].unmatched, 0);
    }

    #[test]
    fn take_completed_record_refuses_a_catalog_that_reuses_an_id() {
        let line = "chr1\t100\t200\tID=TR1;MOTIFS=CAG;STRUC=<TR>\n";
        let other_line = "chr1\t300\t400\tID=TR1;MOTIFS=CAG;STRUC=<TR>\n";
        let mut completed_ids = HashMap::from([("TR1".to_string(), kept_records(1))]);

        // The record kept for TR1 belongs to one of the two loci that share the ID, and nothing in
        // the output says which, so the second line is an error rather than a guess.
        assert!(take_completed_record(&mut completed_ids, line).unwrap());
        let error = take_completed_record(&mut completed_ids, other_line)
            .expect_err("a reused ID that the kept records cannot account for should fail");
        assert!(error.contains("'TR1'"), "{error}");
    }

    #[test]
    fn take_completed_record_keeps_loci_the_previous_run_did_not_reach() {
        let mut completed_ids = HashMap::from([("TR1".to_string(), kept_records(1))]);
        assert!(!take_completed_record(
            &mut completed_ids,
            "chr1\t100\t200\tID=TR2;MOTIFS=CAG;STRUC=<TR>\n"
        )
        .unwrap());
        assert!(!take_completed_record(
            &mut completed_ids,
            "chr1\t100\t200\tMOTIFS=CAG;STRUC=<TR>\n"
        )
        .unwrap());
        assert_eq!(completed_ids["TR1"].unmatched, 1);
    }

    #[test]
    fn check_catalog_ids_accepts_a_catalog_whose_loci_all_have_their_own_id() {
        let dir = tempfile::tempdir().unwrap();
        let catalog = write_catalog(dir.path(), &["TR1", "TR2", "TR3"]);
        let (tr1_start, tr1_end) = catalog_region(0);
        let (tr2_start, tr2_end) = catalog_region(1);
        let mut completed_ids = HashMap::from([
            (
                "TR1".to_string(),
                kept_record_at("chr1", tr1_start, tr1_end),
            ),
            (
                "TR2".to_string(),
                kept_record_at("chr1", tr2_start, tr2_end),
            ),
        ]);

        check_catalog_ids(&catalog, &mut completed_ids).unwrap();

        // The counts are left ready for the run itself to consume them.
        assert_eq!(completed_ids["TR1"].unmatched, 1);
        assert_eq!(completed_ids["TR2"].unmatched, 1);
    }

    #[test]
    fn check_catalog_ids_refuses_a_catalog_with_more_loci_of_an_id_than_there_are_records() {
        let dir = tempfile::tempdir().unwrap();
        let catalog = write_catalog(dir.path(), &["TR1", "TR2", "TR1"]);
        let (start, end) = catalog_region(0);
        let mut completed_ids =
            HashMap::from([("TR1".to_string(), kept_record_at("chr1", start, end))]);

        let error = check_catalog_ids(&catalog, &mut completed_ids)
            .expect_err("a catalog that reuses an ID should fail");
        assert!(error.contains("'TR1'"), "{error}");
        assert!(error.contains("BED line 3"), "{error}");
    }

    #[test]
    fn check_catalog_ids_accepts_a_reused_id_that_every_record_is_there_for() {
        let dir = tempfile::tempdir().unwrap();
        let catalog = write_catalog(dir.path(), &["TR1", "TR1"]);
        let (first_start, first_end) = catalog_region(0);
        let (second_start, second_end) = catalog_region(1);
        let mut kept = kept_record_at("chr1", first_start, first_end);
        kept.add_record(&test_definition("chr1", second_start, second_end));
        let mut completed_ids = HashMap::from([("TR1".to_string(), kept)]);

        check_catalog_ids(&catalog, &mut completed_ids).unwrap();
        assert_eq!(completed_ids["TR1"].unmatched, 2);
    }

    #[test]
    fn check_catalog_ids_refuses_a_catalog_whose_locus_was_given_other_motifs() {
        let dir = tempfile::tempdir().unwrap();
        let catalog = write_catalog(dir.path(), &["TR1", "TR2"]);
        let (start, end) = catalog_region(0);
        let mut kept = KeptRecords::default();
        kept.add_record(&LocusDefinition {
            region: GenomicRegion {
                contig: "chr1".to_string(),
                start,
                end,
            },
            motifs: "CAGG".to_string(),
            struc: "<TR>".to_string(),
        });
        let mut completed_ids = HashMap::from([("TR1".to_string(), kept)]);

        // The locus lies where it did, but it is no longer built of the same motifs, so the
        // genotype the previous run wrote for it is not the genotype this catalog asks for.
        let error = check_catalog_ids(&catalog, &mut completed_ids)
            .expect_err("a locus given other motifs since the previous run should fail");
        assert!(error.contains("'TR1'"), "{error}");
    }

    #[test]
    fn check_catalog_ids_refuses_a_catalog_whose_locus_has_moved() {
        let dir = tempfile::tempdir().unwrap();
        let catalog = write_catalog(dir.path(), &["TR1", "TR2"]);
        let mut completed_ids =
            HashMap::from([("TR1".to_string(), kept_record_at("chr1", 9000, 9030))]);

        let error = check_catalog_ids(&catalog, &mut completed_ids)
            .expect_err("a locus that moved since the previous run should fail");
        assert!(error.contains("'TR1'"), "{error}");
    }

    #[test]
    fn check_catalog_ids_refuses_a_catalog_that_no_longer_holds_the_locus() {
        let dir = tempfile::tempdir().unwrap();
        let catalog = write_catalog(dir.path(), &["TR2", "TR3"]);
        let (start, end) = catalog_region(0);
        let mut completed_ids =
            HashMap::from([("TR1".to_string(), kept_record_at("chr1", start, end))]);

        let error = check_catalog_ids(&catalog, &mut completed_ids)
            .expect_err("a locus taken out of the catalog should fail");
        assert!(error.contains("'TR1'"), "{error}");
    }

    #[test]
    fn parse_locus_id_reads_the_id_field_of_a_catalog_line() {
        assert_eq!(
            parse_locus_id("chr1\t100\t200\tID=TR1;MOTIFS=CAG;STRUC=<TR>\n"),
            Some("TR1")
        );
        assert_eq!(
            parse_locus_id("chr1 100 200 MOTIFS=CAG;ID=TR1;STRUC=<TR>"),
            Some("TR1")
        );
    }

    #[test]
    fn parse_locus_id_returns_nothing_without_an_id_field() {
        assert_eq!(parse_locus_id("chr1\t100\t200\tMOTIFS=CAG"), None);
        assert_eq!(parse_locus_id("chr1\t100\t200\tGENEID=HTT"), None);
        assert_eq!(parse_locus_id("chr1\t100\t200"), None);
        assert_eq!(parse_locus_id(""), None);
    }

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
