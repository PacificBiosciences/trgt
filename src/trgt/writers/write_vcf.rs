//! Defines the `VcfWriter` struct and associated functions for creating and writing results to a VCF file.

use crate::utils::Result;
use crate::{
    trgt::{
        locus::Locus,
        workflows::{Genotype, LocusResult},
    },
    utils::OutputType,
};
use rust_htslib::{
    bam,
    bcf::{
        self,
        record::{GenotypeAllele, Numeric},
        Record,
    },
};
use std::{env, io::Write};

#[cfg(test)]
use rust_htslib::htslib;

const MISSING_FLOAT: f32 = f32::from_bits(0x7F80_0001);

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
    ) -> Result<VcfWriter> {
        let vcf_header = Self::create_header(sample_name, bam_header);

        let output_type = output_type.unwrap_or_default();
        log::debug!("Creating {output_type} writer");
        let (is_uncompressed, format) = output_type.writer_options();
        let writer = bcf::Writer::from_path(output_path, &vcf_header, is_uncompressed, format)
            .map_err(|_| format!("Invalid VCF output path: {output_path}"))?;

        let record = writer.empty_record();

        Ok(VcfWriter { writer, record })
    }

    fn create_header(sample_name: &str, bam_header: &bam::Header) -> bcf::header::Header {
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

        let mut writer =
            VcfWriter::new(tmp.path().to_str().unwrap(), "sample", &bam_header, None).unwrap();

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
}
