//! Defines the `VcfWriter` struct and associated functions for creating and writing results to a VCF file.
//!

use crate::trgt::{
    locus::Locus,
    workflows::{Genotype, LocusResult},
};
use crate::utils::Result;
use itertools::Itertools;
use rust_htslib::{
    bam,
    bcf::{self, record::GenotypeAllele, Format, Record},
};
use std::env;

/// Header lines defining the INFO and FORMAT fields for the VCF file.
const VCF_LINES: [&str; 12] = [
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
];

/// Structure for writing VCF records from genotyping results.
pub struct VcfWriter {
    /// The VCF writer used to write records to the VCF file.
    writer: bcf::Writer,
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
    ) -> Result<VcfWriter> {
        let mut vcf_header = bcf::header::Header::new();

        for line in VCF_LINES.iter() {
            vcf_header.push_record(line.as_bytes());
        }

        if let Some(records) = bam_header.to_hashmap().get("SQ") {
            for record in records {
                let contig_line =
                    format!(r#"##contig=<ID={},length={}>"#, record["SN"], record["LN"]);
                vcf_header.push_record(contig_line.as_bytes());
            }
        }

        let line = format!(
            "##{}Version={}",
            env!("CARGO_PKG_NAME"),
            *crate::cli::FULL_VERSION
        );
        vcf_header.push_record(line.as_bytes());

        let args: Vec<String> = env::args().collect();
        let command_line = args.join(" ");
        let line = format!("##{}Command={}", env!("CARGO_PKG_NAME"), command_line);
        vcf_header.push_record(line.as_bytes());

        vcf_header.push_sample(sample_name.as_bytes());

        let writer = bcf::Writer::from_path(output_path, &vcf_header, false, Format::Vcf)
            .map_err(|_| format!("Invalid VCF output path: {}", output_path))?;

        Ok(VcfWriter { writer })
    }

    /// Writes a VCF record for a given locus and its genotyping results.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `results` - `LocusResult` struct containing genotyping results.
    pub fn write(&mut self, locus: &Locus, results: &LocusResult) {
        let mut record = self.writer.empty_record();
        self.add_locus_info(locus, &mut record);
        if !results.genotype.is_empty() {
            self.add_allele_info(locus, results, &mut record);
        } else {
            self.add_missing_allele_info(locus, &mut record);
        }
        self.writer.write(&record).unwrap();
    }

    /// Adds basic locus information to a VCF record.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `record` - Mutable `Record` struct where the information will be added.
    fn add_locus_info(&mut self, locus: &Locus, record: &mut Record) {
        let contig = locus.region.contig.as_bytes();
        let rid = self.writer.header().name2rid(contig).unwrap();
        record.set_rid(Some(rid));
        record.set_pos(locus.region.start.saturating_sub(1) as i64);

        let id = locus.id.as_bytes();
        record.push_info_string(b"TRID", &[id]).unwrap();
        record
            .push_info_integer(b"END", &[locus.region.end as i32])
            .unwrap();
        let motifs = locus.motifs.join(",");
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
    fn add_missing_allele_info(&mut self, locus: &Locus, record: &mut Record) {
        let pad_base = *locus
            .left_flank
            .as_bytes()
            .last()
            .expect("Empty flanks are not allowed");
        let mut tr_seq = vec![pad_base];
        tr_seq.extend(locus.tr.as_bytes());
        record
            .set_alleles(&[&tr_seq])
            .expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::UnphasedMissing])
            .unwrap();

        record.push_format_string(b"AL", &[".".as_bytes()]).unwrap();
        record
            .push_format_string(b"ALLR", &[".".as_bytes()])
            .unwrap();
        record.push_format_string(b"SD", &[".".as_bytes()]).unwrap();
        record.push_format_string(b"MC", &[".".as_bytes()]).unwrap();
        record.push_format_string(b"MS", &[".".as_bytes()]).unwrap();
        record.push_format_string(b"AP", &[".".as_bytes()]).unwrap();
        record.push_format_string(b"AM", &[".".as_bytes()]).unwrap();
    }

    /// Adds allele information from genotyping results to a VCF record.
    ///
    /// # Arguments
    /// * `locus` - `Locus` struct containing locus information.
    /// * `results` - `LocusResult` struct containing genotyping results.
    /// * `record` - Mutable `Record` struct where the information will be added.
    ///
    /// This method encodes the genotype information and adds it to the VCF record.
    fn add_allele_info(&mut self, locus: &Locus, results: &LocusResult, record: &mut Record) {
        VcfWriter::set_gt(locus, results, record);

        let data = VcfWriter::encode_al(&results.genotype);
        record
            .push_format_string(b"AL", &[data.as_bytes()])
            .unwrap();

        let data = VcfWriter::encode_allr(&results.genotype);
        record
            .push_format_string(b"ALLR", &[data.as_bytes()])
            .unwrap();

        let data = VcfWriter::encode_hd(&results.genotype);
        record
            .push_format_string(b"SD", &[data.as_bytes()])
            .unwrap();

        let data = VcfWriter::encode_mc(&results.genotype);
        record
            .push_format_string(b"MC", &[data.as_bytes()])
            .unwrap();

        let data = VcfWriter::encode_ms(&results.genotype);
        record
            .push_format_string(b"MS", &[data.as_bytes()])
            .unwrap();

        let data = VcfWriter::encode_ap(&results.genotype);
        record
            .push_format_string(b"AP", &[data.as_bytes()])
            .unwrap();

        let data = VcfWriter::encode_am(&results.genotype);
        record
            .push_format_string(b"AM", &[data.as_bytes()])
            .unwrap();
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
    fn set_gt(locus: &Locus, results: &LocusResult, record: &mut Record) {
        let mut seqs = vec![locus.tr.as_bytes()];
        let mut indexes = Vec::new();
        let gt = &results.genotype;
        for allele in gt {
            if allele.seq == locus.tr {
                indexes.push(GenotypeAllele::Unphased(0));
                continue;
            }

            if seqs.len() == 1 {
                indexes.push(GenotypeAllele::Unphased(1));
                seqs.push(allele.seq.as_bytes());
            } else if gt[0].seq == gt[1].seq {
                indexes.push(GenotypeAllele::Unphased(1));
            } else {
                indexes.push(GenotypeAllele::Unphased(2));
                seqs.push(allele.seq.as_bytes());
            }
        }

        let pad_base = *locus
            .left_flank
            .as_bytes()
            .last()
            .expect("Empty flanks are not allowed");
        let padded_seqs = seqs
            .iter()
            .map(|s| {
                let mut padded_seq = vec![pad_base];
                padded_seq.extend(s.iter());
                padded_seq
            })
            .collect_vec();
        let encoding = padded_seqs.iter().map(|s| &s[..]).collect_vec();
        record
            .set_alleles(&encoding)
            .expect("Failed to set alleles");

        record.push_genotypes(&indexes).unwrap();
    }

    /// Encodes the allele lengths as a comma-separated string.
    ///
    /// # Arguments
    /// * `results` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded allele lengths.
    fn encode_al(results: &Genotype) -> String {
        let mut encoding = String::new();
        for hap in results {
            if !encoding.is_empty() {
                encoding += ",";
            }
            encoding += &hap.seq.len().to_string();
        }
        encoding
    }

    /// Encodes the motif counts for each allele.
    ///
    /// # Arguments
    /// * `genotype` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded motif counts.
    fn encode_mc(genotype: &Genotype) -> String {
        let mut encoding = Vec::new();
        for allele in genotype {
            encoding.push(
                allele
                    .annotation
                    .motif_counts
                    .iter()
                    .map(|c| c.to_string())
                    .join("_"),
            );
        }
        encoding.join(",")
    }
    /// Encodes the motif spans for each allele.
    ///
    /// # Arguments
    /// * `results` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded motif spans.
    fn encode_ms(results: &Genotype) -> String {
        let mut encoding = String::new();
        for hap in results {
            let sizes = match &hap.annotation.labels {
                None => vec![".".to_string()],
                Some(value) => value
                    .iter()
                    .map(|s| format!("{}({}-{})", s.motif_index, s.start, s.end))
                    .collect_vec(),
            };
            if !encoding.is_empty() {
                encoding += ",";
            }
            encoding += &sizes.join("_");
        }
        encoding
    }

    /// Encodes the allele purity for each allele.
    ///
    /// # Arguments
    /// * `genotype` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded allele purities.
    fn encode_ap(genotype: &Genotype) -> String {
        genotype
            .iter()
            .map(|a| {
                if a.annotation.purity.is_nan() {
                    ".".to_string()
                } else {
                    format!("{:.6}", a.annotation.purity)
                }
            })
            .join(",")
    }

    /// Encodes the length range for each allele.
    ///
    /// # Arguments
    /// * `results` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded length ranges.
    fn encode_allr(results: &Genotype) -> String {
        let mut encoding = Vec::new();
        for hap in results {
            encoding.push(format!("{}-{}", hap.ci.0, hap.ci.1));
        }
        encoding.join(",")
    }

    /// Encodes the number of spanning reads supporting each allele.
    ///
    /// # Arguments
    /// * `results` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded number of spanning reads.
    fn encode_hd(results: &Genotype) -> String {
        let mut encoding = String::new();
        for hap in results {
            if !encoding.is_empty() {
                encoding += ",";
            }
            encoding += &hap.num_spanning.to_string();
        }
        encoding
    }

    /// Encodes the mean methylation level for each allele.
    ///
    /// # Arguments
    /// * `results` - `Genotype` struct containing the alleles.
    ///
    /// # Returns
    /// Returns a `String` with the encoded mean methylation levels.
    fn encode_am(results: &Genotype) -> String {
        let mut encoding = String::new();
        for hap in results {
            if !encoding.is_empty() {
                encoding += ",";
            }
            encoding += &match hap.meth {
                Some(value) => format!("{:.2}", value),
                None => ".".to_string(),
            };
        }
        encoding
    }
}
