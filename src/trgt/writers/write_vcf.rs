//! Defines the `VcfWriter` struct and associated functions for creating and writing results to a VCF file.

use crate::utils::Result;
use crate::{
    trgt::{
        locus::{check_catalog_ids, KeptRecords, Locus, LocusDefinition},
        workflows::{Genotype, LocusResult},
    },
    utils::{GenomicRegion, InputSource, OutputType},
};
use rust_htslib::{
    bam,
    bcf::{
        self,
        header::HeaderView,
        record::{GenotypeAllele, Numeric},
        Read as BcfRead, Record,
    },
};
use std::{
    collections::{HashMap, HashSet},
    env, fs,
    io::Write,
    path::Path,
};

#[cfg(test)]
use rust_htslib::htslib;

const MISSING_FLOAT: f32 = f32::from_bits(0x7F80_0001);

/// Suffix of the temporary file that a resumed run writes the kept records to before it replaces
/// the output file.
const RESUME_TMP_SUFFIX: &str = ".trgt-resume-tmp";

/// Header lines defining the INFO and FORMAT fields for the VCF file.
const VCF_LINES: [&str; 13] = [
    r#"##INFO=<ID=TRID,Number=1,Type=String,Description="Tandem repeat ID">"#,
    r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">"#,
    r#"##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs that the tandem repeat is composed of">"#,
    r#"##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure of the region">"#,
    r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
    r#"##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">"#,
    r#"##FORMAT=<ID=ALLR,Number=.,Type=String,Description="Length range per allele">"#,
    r#"##FORMAT=<ID=SD,Number=.,Type=Integer,Description="Number of spanning reads supporting per allele">"#,
    r#"##FORMAT=<ID=MC,Number=.,Type=String,Description="Motif counts per allele">"#,
    r#"##FORMAT=<ID=MS,Number=.,Type=String,Description="Motif spans per allele">"#,
    r#"##FORMAT=<ID=AP,Number=.,Type=Float,Description="Allele purity per allele">"#,
    r#"##FORMAT=<ID=AM,Number=.,Type=Float,Description="Mean methylation level per allele">"#,
    r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#,
];

/// Structure for writing VCF records from genotyping results.
pub struct VcfWriter {
    /// The VCF writer used to write records to the VCF file.
    writer: bcf::Writer,
    /// Reusable record to avoid allocation on every write.
    record: Record,
}

impl VcfWriter {
    /// Constructs a new `VcfWriter` instance.
    ///
    /// # Arguments
    /// * `output_path` - Path of the output VCF file.
    /// * `sample_name` - The name of the sample to be written to the VCF file.
    /// * `bam_header` - The BAM header used to extract contig information for the VCF header.
    ///
    /// # Returns
    /// Returns a `Result` with either a new `VcfWriter` instance or an error message.
    pub fn new(
        output_path: &str,
        sample_name: &str,
        bam_header: &bam::Header,
        output_type: Option<OutputType>,
        params: &str,
    ) -> Result<VcfWriter> {
        let vcf_header = Self::create_header(sample_name, bam_header, params);
        Self::open(output_path, &vcf_header, output_type)
    }

    /// Constructs a new `VcfWriter` instance that continues the output VCF of a previous run.
    ///
    /// The records already in `output_path` are copied into a temporary file that then replaces
    /// `output_path`, so that the returned writer adds its records after the ones that were kept.
    /// Copying stops at the first record that cannot be read, which is how a run that was killed
    /// while writing gets picked up. An output file whose header cannot be read at all is discarded
    /// and genotyping starts over from the beginning of the catalog.
    ///
    /// # Arguments
    /// * `output_path` - Path of the output VCF file of the previous run.
    /// * `sample_name` - The name of the sample to be written to the VCF file.
    /// * `bam_header` - The BAM header used to extract contig information for the VCF header.
    /// * `repeats_src` - The repeat catalog, whose IDs the kept records are checked against before
    ///   the output file is replaced.
    ///
    /// # Returns
    /// Returns a `Result` with either the writer and the number of records the previous run left
    /// for each locus ID, or an error message.
    pub fn resume(
        output_path: &str,
        sample_name: &str,
        bam_header: &bam::Header,
        output_type: Option<OutputType>,
        repeats_src: &InputSource,
        params: &str,
    ) -> Result<(VcfWriter, HashMap<String, KeptRecords>)> {
        let mut vcf_header = Self::create_header(sample_name, bam_header, params);

        let Some(mut reader) = open_previous_output(output_path) else {
            return Ok((
                Self::open(output_path, &vcf_header, output_type)?,
                HashMap::new(),
            ));
        };
        check_resume_header(reader.header(), sample_name, bam_header, params)?;

        // The rebuilt header describes the resuming run, so the command lines of the runs that
        // produced the kept records are carried over to keep their provenance in the file.
        for command_line in earlier_command_lines(reader.header()) {
            vcf_header.push_record(command_line.as_bytes());
        }

        let tmp_path = format!("{output_path}{RESUME_TMP_SUFFIX}");
        let mut writer = Self::open(&tmp_path, &vcf_header, output_type)?;
        // The catalog is checked against the kept records before the output file is replaced, so
        // that a catalog the kept records cannot be matched to leaves that file as it was.
        let completed_ids = match copy_and_check_records(&mut reader, &mut writer, repeats_src) {
            Ok(completed_ids) => completed_ids,
            Err(e) => {
                drop(writer);
                let _ = fs::remove_file(&tmp_path);
                return Err(e);
            }
        };

        fs::rename(&tmp_path, output_path).map_err(|e| {
            let _ = fs::remove_file(&tmp_path);
            format!("Failed to replace {output_path} with the resumed output: {e}")
        })?;

        // Logged at warning level because the default verbosity hides anything below it, and a run
        // that asked to resume has to be able to see whether it did.
        log::warn!(
            "Resuming from {output_path} with {} record(s) of the previous run kept",
            completed_ids.values().map(|kept| kept.total).sum::<u32>()
        );
        Ok((writer, completed_ids))
    }

    fn open(
        output_path: &str,
        vcf_header: &bcf::header::Header,
        output_type: Option<OutputType>,
    ) -> Result<VcfWriter> {
        let output_type = output_type.unwrap_or_default();
        log::debug!("Creating {output_type} writer");
        let (is_uncompressed, format) = output_type.writer_options();
        let writer = bcf::Writer::from_path(output_path, vcf_header, is_uncompressed, format)
            .map_err(|_| format!("Invalid VCF output path: {output_path}"))?;

        let record = writer.empty_record();

        Ok(VcfWriter { writer, record })
    }

    fn create_header(
        sample_name: &str,
        bam_header: &bam::Header,
        params: &str,
    ) -> bcf::header::Header {
        let mut vcf_header = bcf::header::Header::new();

        for line in VCF_LINES {
            vcf_header.push_record(line.as_bytes());
        }

        if let Some(records) = bam_header.to_hashmap().get("SQ") {
            for record in records {
                let contig_line = {
                    let name = &record["SN"];
                    let length = &record["LN"];
                    format!(r#"##contig=<ID={name},length={length}>"#)
                };
                vcf_header.push_record(contig_line.as_bytes());
            }
        }

        let package_name = env!("CARGO_PKG_NAME");
        let full_version = crate::cli::FULL_VERSION;
        let version_line = format!("##{package_name}Version={full_version}");
        vcf_header.push_record(version_line.as_bytes());

        let params_line = format!("##{package_name}Params={params}");
        vcf_header.push_record(params_line.as_bytes());

        let command_line = {
            let args: Vec<String> = env::args().collect();
            args.join(" ")
        };
        let command_line = format!("##{package_name}Command={command_line}");
        vcf_header.push_record(command_line.as_bytes());

        vcf_header.push_sample(sample_name.as_bytes());
        vcf_header
    }

    /// Writes a VCF record for a given locus and its genotyping results.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `results` - `LocusResult` struct containing genotyping results.
    pub fn write(&mut self, locus: &Locus, results: &LocusResult) {
        self.record.clear();
        VcfWriter::add_locus_info(&self.writer, locus, &mut self.record);
        if results.genotype.is_empty() {
            VcfWriter::add_missing_allele_info(locus, &mut self.record);
        } else {
            VcfWriter::add_allele_info(locus, results, &mut self.record);
        }
        self.writer.write(&self.record).unwrap();
    }

    /// Adds basic locus information to a VCF record.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `record` - Mutable `Record` struct where the information will be added.
    fn add_locus_info(writer: &bcf::Writer, locus: &Locus, record: &mut Record) {
        let contig = locus.region.contig.as_bytes();
        let rid = writer.header().name2rid(contig).unwrap_or_else(|_| {
            panic!("Contig {:?} not present in VCF header", locus.region.contig)
        });
        record.set_rid(Some(rid));
        record.set_pos(locus.region.start.saturating_sub(1) as i64);
        record.set_qual(MISSING_FLOAT);

        let id = locus.id.as_bytes();
        record.push_info_string(b"TRID", &[id]).unwrap();
        record
            .push_info_integer(b"END", &[locus.region.end as i32])
            .unwrap();
        let motifs = locus
            .motifs
            .iter()
            .map(|m| std::str::from_utf8(m).unwrap())
            .collect::<Vec<_>>()
            .join(",");
        record
            .push_info_string(b"MOTIFS", &[motifs.as_bytes()])
            .unwrap();
        record
            .push_info_string(b"STRUC", &[locus.struc.as_bytes()])
            .unwrap();
    }

    /// Adds allele information for missing genotypes to a VCF record.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `record` - Mutable `Record` struct where the information will be added.
    fn add_missing_allele_info(locus: &Locus, record: &mut Record) {
        let pad_base = *locus
            .left_flank
            .last()
            .expect("Empty flanks are not allowed");
        let mut tr_seq = vec![pad_base];
        tr_seq.extend(&locus.tr);
        record
            .set_alleles(&[&tr_seq])
            .expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::UnphasedMissing])
            .unwrap();

        record
            .push_format_integer(b"AL", &[<i32 as Numeric>::missing()])
            .unwrap();
        record
            .push_format_string(b"ALLR", &[".".as_bytes()])
            .unwrap();
        record
            .push_format_integer(b"SD", &[<i32 as Numeric>::missing()])
            .unwrap();
        record.push_format_string(b"MC", &[".".as_bytes()]).unwrap();
        record.push_format_string(b"MS", &[".".as_bytes()]).unwrap();
        record
            .push_format_float(b"AP", &[<f32 as Numeric>::missing()])
            .unwrap();
        record
            .push_format_float(b"AM", &[<f32 as Numeric>::missing()])
            .unwrap();
        // The PS field is optional; omit it for missing genotypes.
    }

    /// Adds allele information from genotyping results to a VCF record.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `results` - `LocusResult` struct containing genotyping results.
    /// * `record` - Mutable `Record` struct where the information will be added.
    ///
    /// This method encodes the genotype information and adds it to the VCF record.
    fn add_allele_info(locus: &Locus, results: &LocusResult, record: &mut Record) {
        let ps_to_write = VcfWriter::set_gt(locus, results, record);
        record
            .push_format_integer(b"AL", &VcfWriter::encode_al(&results.genotype))
            .unwrap();
        record
            .push_format_string(b"ALLR", &[VcfWriter::encode_allr(&results.genotype)])
            .unwrap();
        record
            .push_format_integer(b"SD", &VcfWriter::encode_hd(&results.genotype))
            .unwrap();
        record
            .push_format_string(b"MC", &[VcfWriter::encode_mc(&results.genotype)])
            .unwrap();
        record
            .push_format_string(b"MS", &[VcfWriter::encode_ms(&results.genotype)])
            .unwrap();
        record
            .push_format_float(b"AP", &VcfWriter::encode_ap(&results.genotype))
            .unwrap();
        record
            .push_format_float(b"AM", &VcfWriter::encode_am(&results.genotype))
            .unwrap();

        if let Some(ps) = ps_to_write {
            record
                .push_format_integer(b"PS", &[ps])
                .expect("Failed to set PS");
        }
    }

    /// This method encodes the genotype information and adds it to the VCF record.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing information about the tandem repeat.
    /// * `results` - `LocusResult` struct containing genotyping results for the locus.
    /// * `record` - A mutable `Record` struct where the genotype information will be added.
    ///
    /// This method constructs the genotype (GT) field for the VCF record by comparing the genotyped alleles
    /// to the reference tandem repeat sequence. It assigns allele indexes and encodes them in VCF format.
    /// It returns the phase set identifier to write to the PS field, if any. The PS field is written after all other format fields.
    fn set_gt(locus: &Locus, results: &LocusResult, record: &mut Record) -> Option<i32> {
        let gt = &results.genotype;
        let mut seqs = vec![locus.tr.as_slice()];
        let mut allele_to_vcf_digit = Vec::with_capacity(gt.len());

        for allele in gt.iter() {
            if allele.seq == locus.tr {
                allele_to_vcf_digit.push(0);
            } else if seqs.len() == 1 {
                seqs.push(allele.seq.as_slice());
                allele_to_vcf_digit.push(1);
            } else if gt.len() == 2 && gt[0].seq == gt[1].seq {
                allele_to_vcf_digit.push(1);
            } else {
                if seqs.len() == 2 {
                    seqs.push(allele.seq.as_slice());
                }
                allele_to_vcf_digit.push(2);
            }
        }

        let pad_base = *locus
            .left_flank
            .last()
            .expect("Empty flanks are not allowed");
        let padded_seqs: Vec<Vec<u8>> = seqs
            .iter()
            .map(|s| {
                let mut v = Vec::with_capacity(s.len() + 1);
                v.push(pad_base);
                v.extend_from_slice(s);
                v
            })
            .collect();
        let encoding: Vec<&[u8]> = padded_seqs.iter().map(Vec::as_slice).collect();
        record
            .set_alleles(&encoding)
            .expect("Failed to set alleles");

        let phase_info = results.phase_result.as_ref().and_then(|pr| {
            if gt.len() != 2 {
                return None;
            }
            let (idx_a, idx_b) = pr.phased_gt_digits;
            let (Some(&a), Some(&b)) = (
                allele_to_vcf_digit.get(idx_a as usize),
                allele_to_vcf_digit.get(idx_b as usize),
            ) else {
                return None;
            };
            (a != b).then_some((a, b, pr.ps))
        });

        let mut gt_field: Vec<GenotypeAllele> = Vec::with_capacity(allele_to_vcf_digit.len());
        let ps_to_write = if let Some((a, b, ps)) = phase_info {
            gt_field.push(GenotypeAllele::Phased(a));
            gt_field.push(GenotypeAllele::Phased(b));
            ps.map(|v| v as i32)
        } else {
            gt_field.extend(
                allele_to_vcf_digit
                    .iter()
                    .copied()
                    .map(GenotypeAllele::Unphased),
            );
            None
        };

        record.push_genotypes(&gt_field).unwrap();
        ps_to_write
    }

    fn encode_al(genotype: &Genotype) -> Vec<i32> {
        genotype
            .iter()
            .map(|allele| {
                i32::try_from(allele.seq.len())
                    .expect("allele length should fit in i32 format field")
            })
            .collect()
    }

    fn encode_allr(genotype: &Genotype) -> Vec<u8> {
        let mut buf = Vec::new();
        for (idx, allele) in genotype.iter().enumerate() {
            if idx > 0 {
                buf.push(b',');
            }
            write!(&mut buf, "{}-{}", allele.ci.0, allele.ci.1).unwrap();
        }
        buf
    }

    fn encode_hd(genotype: &Genotype) -> Vec<i32> {
        genotype
            .iter()
            .map(|allele| {
                i32::try_from(allele.num_spanning)
                    .expect("spanning read count should fit in i32 format field")
            })
            .collect()
    }

    fn encode_mc(genotype: &Genotype) -> Vec<u8> {
        let mut buf = Vec::new();
        for (idx, allele) in genotype.iter().enumerate() {
            if idx > 0 {
                buf.push(b',');
            }
            for (mc_idx, count) in allele.annotation.motif_counts.iter().enumerate() {
                if mc_idx > 0 {
                    buf.push(b'_');
                }
                write!(&mut buf, "{}", count).unwrap();
            }
        }
        buf
    }

    fn encode_ms(genotype: &Genotype) -> Vec<u8> {
        let mut buf = Vec::new();
        for (idx, allele) in genotype.iter().enumerate() {
            if idx > 0 {
                buf.push(b',');
            }
            match &allele.annotation.labels {
                None => buf.push(b'.'),
                Some(v) => {
                    for (lab_idx, s) in v.iter().enumerate() {
                        if lab_idx > 0 {
                            buf.push(b'_');
                        }
                        write!(&mut buf, "{}({}-{})", s.motif_index, s.start, s.end).unwrap();
                    }
                }
            }
        }
        buf
    }

    fn encode_ap(genotype: &Genotype) -> Vec<f32> {
        genotype
            .iter()
            .map(|allele| match allele.annotation.purity {
                x if x.is_nan() => <f32 as Numeric>::missing(),
                x => ((x * 1_000_000.0).round() / 1_000_000.0) as f32,
            })
            .collect()
    }

    fn encode_am(genotype: &Genotype) -> Vec<f32> {
        genotype
            .iter()
            .map(|allele| match allele.meth {
                Some(x) => ((x * 100.0).round() / 100.0) as f32,
                None => <f32 as Numeric>::missing(),
            })
            .collect()
    }

    #[cfg(test)]
    fn record_ptr(&self) -> *mut htslib::bcf1_t {
        self.record.inner
    }
}

/// Opens the output file of a previous run.
///
/// Returns `None` when there is nothing to resume from because the file does not exist, is not a
/// VCF or BCF file, or has no readable header.
fn open_previous_output(output_path: &str) -> Option<bcf::Reader> {
    if !Path::new(output_path).exists() {
        log::warn!("No file at {output_path} to resume from, genotyping the whole catalog");
        return None;
    }

    let reader = match bcf::Reader::from_path(output_path) {
        Ok(reader) => reader,
        Err(e) => {
            log::warn!("Unable to read {output_path} to resume from it ({e}), genotyping the whole catalog and overwriting the file");
            return None;
        }
    };

    if reader.header().inner.is_null() {
        log::warn!("The VCF header of {output_path} could not be read, genotyping the whole catalog and overwriting the file");
        return None;
    }

    Some(reader)
}

/// Copies the records of a previous run into the writer of a resumed run and checks that they can
/// be matched to the loci of the catalog.
fn copy_and_check_records(
    reader: &mut bcf::Reader,
    writer: &mut VcfWriter,
    repeats_src: &InputSource,
) -> Result<HashMap<String, KeptRecords>> {
    let mut completed_ids = copy_kept_records(reader, writer)?;
    check_catalog_ids(repeats_src, &mut completed_ids)?;
    Ok(completed_ids)
}

/// Copies the records of a previous run into the writer of a resumed run.
///
/// Copying stops at the first record that cannot be read or that has no TRID, because the rest of
/// the file cannot be trusted once a record is unreadable.
///
/// The last record of the file is never kept. A run that was killed while writing leaves a record
/// that stops in the middle, and htslib reads a record whose text stops inside the sample column as
/// a complete one: it has the right number of columns, so the fields that the record no longer
/// reaches are simply read as missing. Such a record would be copied with a genotype that is
/// silently wrong. Only the last record of a file can be cut short, so a record is copied once
/// something has been read after it, whether or not that something could be made sense of, and the
/// record left over at the end is dropped. That costs one locus of genotyping per resumed run.
///
/// Returns how many records were kept for each locus ID, so that the catalog reader can skip a
/// locus exactly as many times as the file already holds a record for it.
fn copy_kept_records(
    reader: &mut bcf::Reader,
    writer: &mut VcfWriter,
) -> Result<HashMap<String, KeptRecords>> {
    let mut completed_ids: HashMap<String, KeptRecords> = HashMap::new();
    let mut num_kept = 0;
    let mut pending: Option<(Record, String, LocusDefinition)> = None;
    for result in reader.records() {
        let record = match result {
            Ok(record) => record,
            Err(e) => {
                // Reaching a record at all, even one that cannot be read, proves that the record
                // held back before it is whole, so it is kept.
                keep_pending(&mut pending, writer, &mut completed_ids, &mut num_kept)?;
                log::warn!("Keeping the first {num_kept} record(s) of the previous run, the record after them could not be read ({e}), so its locus and all loci after it are genotyped again");
                break;
            }
        };
        // htslib hands back a record that it could only partly make sense of, for instance one on a
        // contig that the header does not define, with an error code set. Copying such a record
        // would abort the process, so stop at it like at a record that could not be read at all.
        if record.inner().errcode != 0 {
            keep_pending(&mut pending, writer, &mut completed_ids, &mut num_kept)?;
            log::warn!("Keeping the first {num_kept} record(s) of the previous run, the record after them is malformed, so its locus and all loci after it are genotyped again");
            break;
        }

        let Some(id) = record_trid(&record) else {
            keep_pending(&mut pending, writer, &mut completed_ids, &mut num_kept)?;
            log::warn!("Keeping the first {num_kept} record(s) of the previous run, the record after them has no TRID field, so its locus and all loci after it are genotyped again");
            break;
        };
        let Some(definition) = record_definition(&record) else {
            keep_pending(&mut pending, writer, &mut completed_ids, &mut num_kept)?;
            log::warn!("Keeping the first {num_kept} record(s) of the previous run, the locus of the record after them could not be read, so that locus and all loci after it are genotyped again");
            break;
        };

        // Reading this record proved that the record before it is whole.
        keep_pending(&mut pending, writer, &mut completed_ids, &mut num_kept)?;
        pending = Some((record, id, definition));
    }

    let num_duplicated_ids = completed_ids.values().filter(|kept| kept.total > 1).count();
    if num_duplicated_ids > 0 {
        log::warn!("The previous run wrote more than one record for {num_duplicated_ids} of the locus IDs, the catalog may contain duplicate IDs");
    }
    Ok(completed_ids)
}

/// Copies the record held back, if there is one, now that it is known to be whole.
fn keep_pending(
    pending: &mut Option<(Record, String, LocusDefinition)>,
    writer: &mut VcfWriter,
    completed_ids: &mut HashMap<String, KeptRecords>,
    num_kept: &mut usize,
) -> Result<()> {
    let Some((record, id, definition)) = pending.take() else {
        return Ok(());
    };
    copy_record(writer, record, &id)?;
    completed_ids.entry(id).or_default().add_record(&definition);
    *num_kept += 1;
    Ok(())
}

/// Writes a record of a previous run into the output of a resumed run.
fn copy_record(writer: &mut VcfWriter, mut record: Record, id: &str) -> Result<()> {
    writer.writer.translate(&mut record);
    writer
        .writer
        .write(&record)
        .map_err(|e| format!("Failed to copy the record of locus {id} into the output: {e}"))
}

/// Returns what the catalog of the previous run defined a record's locus as, or `None` if the
/// record does not say.
///
/// Every part is read the way the catalog gives it: the contig, the start that TRGT wrote as the
/// record's position, the end that it wrote as the END field, and the motifs and structure it wrote
/// as the MOTIFS and STRUC fields.
fn record_definition(record: &Record) -> Option<LocusDefinition> {
    let contig = record.header().rid2name(record.rid()?).ok()?;
    let end = record.info(b"END").integer().ok()??.first().copied()?;
    Some(LocusDefinition {
        region: GenomicRegion {
            contig: String::from_utf8_lossy(contig).into_owned(),
            start: u32::try_from(record.pos() + 1).ok()?,
            end: u32::try_from(end).ok()?,
        },
        motifs: record_info_string(record, b"MOTIFS")?,
        struc: record_info_string(record, b"STRUC")?,
    })
}

/// Returns the value of one of a record's text INFO fields, spelled the way the catalog spells it.
///
/// htslib hands a value that holds commas back as several strings, one per comma-separated part, so
/// the parts are joined again to give back the value the catalog held. A locus built of several
/// motifs, whose MOTIFS field reads `CAG,CAA`, comes back whole this way.
fn record_info_string(record: &Record, tag: &[u8]) -> Option<String> {
    let value = record.info(tag).string().ok()??;
    let parts: Vec<String> = value
        .iter()
        .map(|part| String::from_utf8_lossy(part).into_owned())
        .collect();
    (!parts.is_empty()).then(|| parts.join(","))
}

/// Returns the TRID of a record, or `None` if the record has no TRID field.
fn record_trid(record: &Record) -> Option<String> {
    record_info_string(record, b"TRID")
}

/// Checks that the output file of a previous run describes the same sample and reference as the
/// run that is resuming it, and that all of its INFO and FORMAT fields are still defined.
///
/// Resuming rewrites the header, so a file from a different sample, reference, or TRGT version
/// would end up with records that no longer match the header that describes them.
fn check_resume_header(
    header: &HeaderView,
    sample_name: &str,
    bam_header: &bam::Header,
    params: &str,
) -> Result<()> {
    let samples = header.samples();
    match samples.as_slice() {
        [sample] if *sample == sample_name.as_bytes() => {}
        [sample] => {
            return Err(format!(
                "Cannot resume: the output file holds genotypes of sample '{}' but this run genotypes sample '{sample_name}'",
                String::from_utf8_lossy(sample)
            ))
        }
        samples => {
            return Err(format!(
                "Cannot resume: the output file holds genotypes of {} samples, TRGT writes one",
                samples.len()
            ))
        }
    }

    let version = header_value(header, &format!("{}Version", env!("CARGO_PKG_NAME")));
    match version.as_deref() {
        Some(version) if version == crate::cli::FULL_VERSION => {}
        Some(version) => {
            return Err(format!(
                "Cannot resume: the output file was written by {package} {version} but this run is {package} {current}, and genotypes of two versions do not belong in one file",
                package = env!("CARGO_PKG_NAME"),
                current = crate::cli::FULL_VERSION
            ))
        }
        None => {
            return Err(format!(
                "Cannot resume: the output file does not say which version of {} wrote it",
                env!("CARGO_PKG_NAME")
            ))
        }
    }

    let previous_params = header_value(header, &format!("{}Params", env!("CARGO_PKG_NAME")));
    match previous_params.as_deref() {
        Some(previous_params) if previous_params == params => {}
        Some(previous_params) => {
            return Err(format!(
                "Cannot resume: the output file was genotyped with settings that differ from this run's.\n  output file: {previous_params}\n  this run:    {params}"
            ))
        }
        None => {
            return Err(format!(
                "Cannot resume: the output file does not say which settings it was genotyped with, so it cannot be told whether they match this run's ({params})"
            ))
        }
    }

    let contig_lengths = bam_contig_lengths(bam_header);
    let defined_ids = defined_field_ids();
    for record in header.header_records() {
        match record {
            bcf::header::HeaderRecord::Contig { values, .. } => {
                let (Some(name), Some(length)) = (values.get("ID"), values.get("length")) else {
                    continue;
                };
                match contig_lengths.get(name.as_str()) {
                    Some(bam_length) if bam_length == length => {}
                    Some(bam_length) => {
                        return Err(format!(
                            "Cannot resume: the output file gives contig {name} a length of {length} but the reads give it a length of {bam_length}"
                        ))
                    }
                    None => {
                        return Err(format!(
                            "Cannot resume: the output file gives a contig {name}, which is missing from the reads"
                        ))
                    }
                }
            }
            bcf::header::HeaderRecord::Info { values, .. }
            | bcf::header::HeaderRecord::Format { values, .. } => {
                let Some(id) = values.get("ID") else { continue };
                if !defined_ids.contains(id.as_str()) {
                    return Err(format!(
                        "Cannot resume: the output file uses the field {id}, which this version of TRGT does not write"
                    ));
                }
            }
            _ => {}
        }
    }

    Ok(())
}

/// Returns the value of a header line such as `##trgtVersion=...`, if the header has one.
fn header_value(header: &HeaderView, key: &str) -> Option<String> {
    header
        .header_records()
        .into_iter()
        .find_map(|record| match record {
            bcf::header::HeaderRecord::Generic {
                key: record_key,
                value,
            } if record_key == key => Some(value),
            _ => None,
        })
}

/// Returns the command lines of the runs that wrote an output file, as header lines.
///
/// The command line of the run that is resuming is left out, since the rebuilt header has it
/// already.
fn earlier_command_lines(header: &HeaderView) -> Vec<String> {
    let key = format!("{}Command", env!("CARGO_PKG_NAME"));
    let current_command = env::args().collect::<Vec<_>>().join(" ");
    let mut command_lines: Vec<String> = Vec::new();
    for record in header.header_records() {
        let bcf::header::HeaderRecord::Generic {
            key: record_key,
            value,
        } = record
        else {
            continue;
        };
        if record_key != key || value == current_command {
            continue;
        }
        let command_line = format!("##{key}={value}");
        if !command_lines.contains(&command_line) {
            command_lines.push(command_line);
        }
    }
    command_lines
}

/// Returns the name and length of every contig in a BAM header.
fn bam_contig_lengths(bam_header: &bam::Header) -> HashMap<String, String> {
    let mut lengths = HashMap::new();
    if let Some(records) = bam_header.to_hashmap().get("SQ") {
        for record in records {
            if let (Some(name), Some(length)) = (record.get("SN"), record.get("LN")) {
                lengths.insert(name.clone(), length.clone());
            }
        }
    }
    lengths
}

/// Returns the IDs of the INFO and FORMAT fields that TRGT defines in its VCF header.
fn defined_field_ids() -> HashSet<&'static str> {
    VCF_LINES
        .iter()
        .filter_map(|line| {
            let values = line.split_once("<ID=")?.1;
            Some(values.split_once(',')?.0)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        trgt::workflows::{Allele, PhaseResult},
        utils::test_util::TestVcfBuilder,
    };
    use rust_htslib::{
        bam::{self, header::HeaderRecord},
        bcf::{self, Read as BcfRead, Reader, Writer},
    };
    use std::str::FromStr;
    use tempfile::NamedTempFile;

    fn empty_record() -> Record {
        let tmp = NamedTempFile::new().unwrap();
        let header = TestVcfBuilder::new()
            .contig("chr1", 10)
            .add_format("PS", "1", "Integer", "Phase set identifier")
            .sample("sample")
            .build_header();
        let writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();
        writer.empty_record()
    }

    #[test]
    fn test_phased_het_writes_phased_gt_and_ps() {
        let locus = Locus::test_base();
        let mut genotype = Genotype::new();
        genotype.push(Allele::test_base(&locus.tr));
        genotype.push(Allele::test_base(b"TT"));

        let results = LocusResult {
            genotype,
            spanning_reads: Vec::new(),
            phase_result: Some(PhaseResult {
                phased_gt_digits: (0, 1),
                ps: Some(10),
            }),
        };

        let mut record = empty_record();
        if let Some(ps) = VcfWriter::set_gt(&locus, &results, &mut record) {
            record
                .push_format_integer(b"PS", &[ps])
                .expect("Failed to set PS");
        }

        let gt = record.genotypes().unwrap();
        assert_eq!("0|1", format!("{}", gt.get(0)));
        let ps = record.format(b"PS").integer().unwrap();
        assert_eq!(ps[0][0], 10);
    }

    #[test]
    fn test_hom_remains_unphased() {
        let locus = Locus::test_base();
        let mut genotype = Genotype::new();
        genotype.push(Allele::test_base(b"TT"));
        genotype.push(Allele::test_base(b"TT"));
        let results = LocusResult {
            genotype,
            spanning_reads: Vec::new(),
            phase_result: Some(PhaseResult {
                phased_gt_digits: (0, 1),
                ps: Some(7),
            }),
        };

        let mut record = empty_record();
        if let Some(ps) = VcfWriter::set_gt(&locus, &results, &mut record) {
            record
                .push_format_integer(b"PS", &[ps])
                .expect("Failed to set PS");
        }

        let gt = record.genotypes().unwrap();
        assert_eq!("1/1", format!("{}", gt.get(0)));
        assert!(record.format(b"PS").integer().is_err());
    }

    #[test]
    fn test_cached_record_reuse() {
        let tmp = NamedTempFile::new().unwrap();
        let mut bam_header = bam::Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 10);
        bam_header.push_record(&sq);

        let mut writer = VcfWriter::new(
            tmp.path().to_str().unwrap(),
            "sample",
            &bam_header,
            None,
            TEST_PARAMS,
        )
        .unwrap();

        let locus = Locus::test_base();
        let mut genotype = Genotype::new();
        genotype.push(Allele::test_base(&locus.tr));
        genotype.push(Allele::test_base(b"TT"));
        let phased = LocusResult {
            genotype,
            spanning_reads: Vec::new(),
            phase_result: Some(PhaseResult {
                phased_gt_digits: (0, 1),
                ps: Some(5),
            }),
        };

        let first_ptr = writer.record_ptr();
        writer.write(&locus, &phased);
        let after_first_ptr = writer.record_ptr();

        let mut second_locus = Locus::test_base();
        second_locus.id = "locus2".into();
        second_locus.region.start = 3;
        second_locus.region.end = 4;
        let mut second_genotype = Genotype::new();
        second_genotype.push(Allele::test_base(&second_locus.tr));
        second_genotype.push(Allele::test_base(b"TTT"));
        let unphased = LocusResult {
            genotype: second_genotype,
            spanning_reads: Vec::new(),
            phase_result: None,
        };

        writer.write(&second_locus, &unphased);
        let after_second_ptr = writer.record_ptr();
        drop(writer);

        let mut reader = Reader::from_path(tmp.path()).unwrap();
        let mut records = reader.records();
        let first = records.next().unwrap().unwrap();
        let second = records.next().unwrap().unwrap();
        assert!(first.format(b"PS").integer().is_ok());
        assert!(second.format(b"PS").integer().is_err());
        assert!(records.next().is_none());

        assert_eq!(first_ptr, after_first_ptr);
        assert_eq!(first_ptr, after_second_ptr);
    }

    /// Stand-in for the genotyping settings a run records in its output.
    const TEST_PARAMS: &str = "genotyper=size;karyotype=XX";

    fn test_bam_header() -> bam::Header {
        let mut bam_header = bam::Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 10);
        bam_header.push_record(&sq);
        bam_header
    }

    /// Uncompressed output keeps the tests that damage the output file readable.
    fn uncompressed_vcf() -> Option<OutputType> {
        Some(OutputType::Vcf {
            is_uncompressed: true,
            level: None,
        })
    }

    fn test_locus(id: &str, start: u32) -> Locus {
        let mut locus = Locus::test_base();
        locus.id = id.into();
        locus.region.start = start;
        locus.region.end = start + 1;
        locus
    }

    fn test_result(locus: &Locus) -> LocusResult {
        let mut genotype = Genotype::new();
        genotype.push(Allele::test_base(&locus.tr));
        genotype.push(Allele::test_base(b"TT"));
        LocusResult {
            genotype,
            spanning_reads: Vec::new(),
            phase_result: None,
        }
    }

    fn write_loci(path: &str, bam_header: &bam::Header, ids: &[&str]) {
        let mut writer =
            VcfWriter::new(path, "sample", bam_header, uncompressed_vcf(), TEST_PARAMS).unwrap();
        for (index, id) in ids.iter().enumerate() {
            let locus = test_locus(id, index as u32 + 1);
            writer.write(&locus, &test_result(&locus));
        }
    }

    /// The kept-record counts of a resume that kept one record for each of `ids`.
    ///
    /// The digests are zero because the catalog the resume checked against describes the same
    /// regions as the kept records, so the two sides cancelled each other out.
    fn kept_once(ids: &[&str]) -> HashMap<String, KeptRecords> {
        ids.iter()
            .map(|id| {
                let kept = KeptRecords {
                    total: 1,
                    unmatched: 1,
                    digest: 0,
                };
                ((*id).to_string(), kept)
            })
            .collect()
    }

    /// Writes a catalog holding one locus per ID, as the source a resumed run checks against.
    fn write_catalog(dir: &Path, ids: &[&str]) -> InputSource {
        let path = dir.join("repeats.bed");
        let lines: String = ids
            .iter()
            .enumerate()
            .map(|(index, id)| {
                let start = index + 1;
                // The motifs and the structure have to be the ones `test_locus` gives its loci,
                // since a resumed run checks those against the records it keeps.
                format!("chr1\t{start}\t{}\tID={id};MOTIFS=T;STRUC=1\n", start + 1)
            })
            .collect();
        fs::write(&path, lines).unwrap();
        InputSource::from_str(path.to_str().unwrap()).unwrap()
    }

    fn read_trids(path: &str) -> Vec<String> {
        let mut reader = Reader::from_path(path).unwrap();
        reader
            .records()
            .map(|record| record_trid(&record.unwrap()).unwrap())
            .collect()
    }

    #[test]
    fn resume_keeps_previous_records_and_adds_new_ones_after_them() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2", "locus3"]);

        let (mut writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        assert_eq!(completed_ids, kept_once(&["locus1", "locus2"]));

        let locus = test_locus("locus3", 3);
        writer.write(&locus, &test_result(&locus));
        drop(writer);

        assert_eq!(read_trids(path), ["locus1", "locus2", "locus3"]);
        assert!(!Path::new(&format!("{path}{RESUME_TMP_SUFFIX}")).exists());
    }

    #[test]
    fn resume_never_keeps_the_last_record_of_the_previous_run() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2"]);

        // The last record is the one a killed run can have left half written, so it is dropped and
        // its locus is genotyped again even when the previous run did finish writing it.
        let (writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        drop(writer);

        assert_eq!(completed_ids, kept_once(&["locus1"]));
        assert_eq!(read_trids(path), ["locus1"]);
    }

    #[test]
    fn resume_drops_a_record_cut_short_inside_its_sample_column() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2", "locus3"]);

        // Cut the last record inside its sample column. htslib reads such a line as a complete
        // record, with the fields it no longer reaches read as missing, so nothing about the record
        // itself gives away that its genotype is only partly written.
        let contents = fs::read_to_string(path).unwrap();
        let last_line_start = contents[..contents.len() - 1].rfind('\n').unwrap() + 1;
        let sample_start = last_line_start + contents[last_line_start..].rfind('\t').unwrap() + 1;
        fs::write(path, &contents[..sample_start + 3]).unwrap();

        let (writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        drop(writer);

        // The half-written record is neither kept in the output nor counted as genotyped, so its
        // locus is genotyped again instead of keeping a genotype that is only partly written.
        assert_eq!(completed_ids, kept_once(&["locus1", "locus2"]));
        assert_eq!(read_trids(path), ["locus1", "locus2"]);
    }

    #[test]
    fn resume_genotypes_everything_when_there_is_no_previous_output() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();

        let (mut writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        assert!(completed_ids.is_empty());

        let locus = test_locus("locus1", 1);
        writer.write(&locus, &test_result(&locus));
        drop(writer);

        assert_eq!(read_trids(path), ["locus1"]);
    }

    #[test]
    fn resume_overwrites_an_output_file_that_is_not_a_vcf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        fs::write(path, b"the previous run wrote nonsense\n").unwrap();

        let (mut writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        assert!(completed_ids.is_empty());

        let locus = test_locus("locus1", 1);
        writer.write(&locus, &test_result(&locus));
        drop(writer);

        assert_eq!(read_trids(path), ["locus1"]);
    }

    #[test]
    fn resume_overwrites_an_output_file_with_a_truncated_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1"]);

        // Keep the first header line only, so that the file still looks like a VCF but has no
        // #CHROM line to end its header.
        let contents = fs::read_to_string(path).unwrap();
        fs::write(path, contents.lines().next().unwrap()).unwrap();

        let (mut writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        assert!(completed_ids.is_empty());

        let locus = test_locus("locus2", 2);
        writer.write(&locus, &test_result(&locus));
        drop(writer);

        assert_eq!(read_trids(path), ["locus2"]);
    }

    #[test]
    fn resume_keeps_the_records_before_a_record_that_cannot_be_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2", "locus3"]);

        // Cut the file before the sample column of the last record, which leaves a line with too
        // few columns for htslib to read.
        let contents = fs::read_to_string(path).unwrap();
        let last_record_start = contents.rfind("chr1\t3\t").unwrap();
        fs::write(path, &contents[..last_record_start + "chr1\t3\t".len()]).unwrap();

        // The record before the unreadable one is kept: bytes follow it, which proves it is whole,
        // so only the locus of the unreadable record has to be genotyped again.
        let (mut writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        assert_eq!(completed_ids, kept_once(&["locus1", "locus2"]));

        let locus = test_locus("locus3", 3);
        writer.write(&locus, &test_result(&locus));
        drop(writer);

        assert_eq!(read_trids(path), ["locus1", "locus2", "locus3"]);
    }

    #[test]
    fn resume_keeps_the_records_before_a_malformed_record() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2"]);

        // Add a record on a contig that the header does not define. htslib reads such a record but
        // flags it, and copying a flagged record would abort the process.
        let contents = fs::read_to_string(path).unwrap();
        let malformed = contents.lines().last().unwrap().replacen("chr1", "chrZ", 1);
        fs::write(path, format!("{contents}{malformed}\n")).unwrap();

        let (mut writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        assert_eq!(completed_ids, kept_once(&["locus1", "locus2"]));

        let locus = test_locus("locus3", 3);
        writer.write(&locus, &test_result(&locus));
        drop(writer);

        assert_eq!(read_trids(path), ["locus1", "locus2", "locus3"]);
    }

    #[test]
    fn resume_refuses_an_output_file_of_another_sample() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1"]);

        let error = VcfWriter::resume(
            path,
            "other_sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .err()
        .expect("resuming the output of another sample should fail");
        assert!(error.contains("sample 'sample'"), "{error}");
        assert_eq!(read_trids(path), ["locus1"]);
    }

    #[test]
    fn resume_refuses_an_output_file_of_another_reference() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        write_loci(path, &test_bam_header(), &["locus1"]);

        let mut other_header = bam::Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 20);
        other_header.push_record(&sq);

        let error = VcfWriter::resume(
            path,
            "sample",
            &other_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .err()
        .expect("resuming the output of another reference should fail");
        assert!(error.contains("length"), "{error}");
    }

    #[test]
    fn resume_keeps_a_locus_built_of_several_motifs() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let bam_header = test_bam_header();

        // A locus whose MOTIFS field holds more than one motif, as most catalogs have. htslib hands
        // such a field back in parts, and a resumed run has to put it together the way the catalog
        // spells it, or it would take an unchanged catalog for a changed one.
        let catalog_path = dir.path().join("multi_motif.bed");
        fs::write(
            &catalog_path,
            "chr1\t1\t2\tID=locus1;MOTIFS=T,TT;STRUC=1\nchr1\t2\t3\tID=locus2;MOTIFS=T,TT;STRUC=1\n",
        )
        .unwrap();
        let catalog = InputSource::from_str(catalog_path.to_str().unwrap()).unwrap();

        let mut writer =
            VcfWriter::new(path, "sample", &bam_header, uncompressed_vcf(), TEST_PARAMS).unwrap();
        for (index, id) in ["locus1", "locus2"].iter().enumerate() {
            let mut locus = test_locus(id, index as u32 + 1);
            locus.motifs = vec![b"T".to_vec(), b"TT".to_vec()];
            writer.write(&locus, &test_result(&locus));
        }
        drop(writer);

        let (writer, completed_ids) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        drop(writer);

        assert_eq!(completed_ids, kept_once(&["locus1"]));
    }

    #[test]
    fn resume_refuses_a_catalog_that_reuses_an_id_without_touching_the_output() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus1", "locus2"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2"]);
        let before = fs::read_to_string(path).unwrap();

        let error = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .err()
        .expect("a catalog that reuses an ID should fail");
        assert!(error.contains("'locus1'"), "{error}");

        // The refusal has to leave the previous run's output as it was, down to the last record.
        assert_eq!(fs::read_to_string(path).unwrap(), before);
        assert!(!Path::new(&format!("{path}{RESUME_TMP_SUFFIX}")).exists());
    }

    #[test]
    fn resume_refuses_an_output_file_of_another_trgt_version() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2"]);

        let contents = fs::read_to_string(path).unwrap();
        let older_version = contents.replace(
            &format!("##trgtVersion={}", crate::cli::FULL_VERSION),
            "##trgtVersion=0.0.0",
        );
        fs::write(path, older_version).unwrap();

        let error = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .err()
        .expect("resuming the output of another version should fail");
        assert!(error.contains("0.0.0"), "{error}");
    }

    #[test]
    fn resume_refuses_an_output_file_genotyped_with_other_settings() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2"]);

        let other_settings = "genotyper=cluster;karyotype=XY";
        let error = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            other_settings,
        )
        .err()
        .expect("resuming an output genotyped with other settings should fail");
        assert!(error.contains(TEST_PARAMS), "{error}");
        assert!(error.contains(other_settings), "{error}");
    }

    #[test]
    fn resume_keeps_the_command_line_of_the_run_that_wrote_the_kept_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1", "locus2"]);

        let earlier_command = "##trgtCommand=trgt genotype --from-an-earlier-run";
        let contents = fs::read_to_string(path).unwrap();
        let current_command = contents
            .lines()
            .find(|line| line.starts_with("##trgtCommand="))
            .unwrap()
            .to_string();
        fs::write(path, contents.replace(&current_command, earlier_command)).unwrap();

        let (writer, _) = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .unwrap();
        drop(writer);

        let resumed = fs::read_to_string(path).unwrap();
        let command_lines: Vec<&str> = resumed
            .lines()
            .filter(|line| line.starts_with("##trgtCommand="))
            .collect();
        assert!(
            command_lines.contains(&earlier_command),
            "{command_lines:?}"
        );
        assert!(
            command_lines.contains(&current_command.as_str()),
            "{command_lines:?}"
        );
    }

    #[test]
    fn resume_refuses_an_output_file_with_an_unknown_field() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.vcf");
        let path = path.to_str().unwrap();
        let catalog = write_catalog(dir.path(), &["locus1", "locus2", "locus3"]);
        let bam_header = test_bam_header();
        write_loci(path, &bam_header, &["locus1"]);

        // Add a field that this version of TRGT does not define, the way an output file written by
        // another version of TRGT would.
        let contents = fs::read_to_string(path).unwrap();
        let (fileformat, rest) = contents.split_once('\n').unwrap();
        let extra_field = r#"##FORMAT=<ID=XX,Number=1,Type=Integer,Description="Dropped field">"#;
        fs::write(path, format!("{fileformat}\n{extra_field}\n{rest}")).unwrap();

        let error = VcfWriter::resume(
            path,
            "sample",
            &bam_header,
            uncompressed_vcf(),
            &catalog,
            TEST_PARAMS,
        )
        .err()
        .expect("resuming an output file with an unknown field should fail");
        assert!(error.contains("XX"), "{error}");
    }

    #[test]
    fn defined_field_ids_holds_every_field_of_the_header() {
        let ids = defined_field_ids();
        assert_eq!(ids.len(), VCF_LINES.len());
        assert!(ids.contains("TRID"));
        assert!(ids.contains("PS"));
    }
}
