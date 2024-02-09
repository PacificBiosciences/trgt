use crate::locus::Locus;
use crate::workflows::{Genotype, LocusResult};
use itertools::Itertools;
use rust_htslib::bam::{self};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self, Format, Record};
use std::env;

pub struct VcfWriter {
    writer: bcf::Writer,
}

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

impl VcfWriter {
    pub fn new(
        output_path: &str,
        sample_name: &str,
        bam_header: &bam::Header,
    ) -> Result<VcfWriter, String> {
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

    pub fn write(&mut self, locus: &Locus, results: &LocusResult) {
        let mut record = self.writer.empty_record();
        self.write_info_fields(locus, results, &mut record);
        if !results.genotype.is_empty() {
            self.write_genotype_fields(locus, results, &mut record);
        } else {
            self.write_missing_genotype_fields(locus, &mut record);
        }
        self.writer.write(&record).unwrap();
    }

    fn write_info_fields(&mut self, locus: &Locus, results: &LocusResult, record: &mut Record) {
        let contig = locus.region.contig.as_bytes();
        let rid = self.writer.header().name2rid(contig).unwrap();
        record.set_rid(Some(rid));

        let no_empty_alleles = results.genotype.iter().all(|a| !a.seq.is_empty());
        record.set_pos(if no_empty_alleles {
            locus.region.start as i64
        } else {
            locus.region.start as i64 - 1
        });

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

    fn write_missing_genotype_fields(&mut self, locus: &Locus, record: &mut Record) {
        record
            .set_alleles(&[locus.tr.as_bytes()])
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

    fn write_genotype_fields(&mut self, locus: &Locus, results: &LocusResult, record: &mut Record) {
        set_gt(locus, results, record);

        let data = encode_al(&results.genotype);
        record
            .push_format_string(b"AL", &[data.as_bytes()])
            .unwrap();

        let data = encode_allr(&results.genotype);
        record
            .push_format_string(b"ALLR", &[data.as_bytes()])
            .unwrap();

        let data = encode_hd(&results.genotype);
        record
            .push_format_string(b"SD", &[data.as_bytes()])
            .unwrap();

        let data = encode_mc(&results.genotype);
        record
            .push_format_string(b"MC", &[data.as_bytes()])
            .unwrap();

        let data = encode_ms(&results.genotype);
        record
            .push_format_string(b"MS", &[data.as_bytes()])
            .unwrap();

        let data = encode_ap(&results.genotype);
        record
            .push_format_string(b"AP", &[data.as_bytes()])
            .unwrap();

        let data = encode_am(&results.genotype);
        record
            .push_format_string(b"AM", &[data.as_bytes()])
            .unwrap();
    }
}

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

    let no_empty_alleles = results.genotype.iter().all(|a| !a.seq.is_empty());
    if no_empty_alleles {
        record.set_alleles(&seqs).expect("Failed to set alleles");
    } else {
        let pad_base = *locus.left_flank.as_bytes().last().unwrap();
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
    }

    record.push_genotypes(&indexes).unwrap();
}

fn encode_al(diplotype: &Genotype) -> String {
    let mut encoding = String::new();

    for hap in diplotype {
        if !encoding.is_empty() {
            encoding += ",";
        }
        encoding += &hap.seq.len().to_string();
    }

    encoding
}

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

fn encode_ms(diplotype: &Genotype) -> String {
    let mut encoding = String::new();
    for hap in diplotype {
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

fn encode_allr(diplotype: &Genotype) -> String {
    let mut encoding = Vec::new();
    for hap in diplotype {
        encoding.push(format!("{}-{}", hap.ci.0, hap.ci.1));
    }

    encoding.join(",")
}

fn encode_hd(diplotype: &Genotype) -> String {
    let mut encoding = String::new();
    for hap in diplotype {
        if !encoding.is_empty() {
            encoding += ",";
        }
        encoding += &hap.num_spanning.to_string();
    }
    encoding
}

fn encode_am(diplotype: &Genotype) -> String {
    let mut encoding = String::new();
    for hap in diplotype {
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
