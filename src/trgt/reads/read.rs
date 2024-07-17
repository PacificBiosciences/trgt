//! Module for representing and building read information from alignment records.
//!

use super::{cigar::Cigar, meth, snp};
use crate::utils::GenomicRegion;
use itertools::Itertools;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::str;

/// Represents methylation information extracted from a read.
///
/// # Attributes
/// * `poses` - Vector of positions where methylation is detected.
/// * `probs` - Vector of probabilities associated with the methylation calls.
#[derive(Debug)]
pub struct MethInfo {
    pub poses: Vec<usize>,
    pub probs: Vec<u8>,
}

/// Represents a single HiFi read from an alignment record.
#[derive(PartialEq, Clone)]
pub struct HiFiRead {
    /// Unique identifier for the read.
    pub id: String,
    /// Flag indicating if the read is from the reverse strand.
    pub is_reverse: bool,
    /// Vector of bases (nucleotides) in the read.
    pub bases: Vec<u8>,
    /// Vector of quality scores for the bases.
    pub quals: Vec<u8>,
    /// Optional vector of methylation calls.
    pub meth: Option<Vec<u8>>,
    /// Optional overall quality score for the read.
    pub read_qual: Option<f64>,
    /// Optional vector of offsets where mismatches occur relative to a locus region.
    pub mismatch_offsets: Option<Vec<i32>>,
    /// Offset from the start of the reference region.
    pub start_offset: i32,
    /// Offset from the end of the reference region.
    pub end_offset: i32,
    /// Optional CIGAR string representing the alignment.
    pub cigar: Option<Cigar>,
    /// Optional haplotype tag.
    pub hp_tag: Option<u8>,
    /// Mapping quality score.
    pub mapq: u8,
}

impl std::fmt::Debug for HiFiRead {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let meth = match &self.meth {
            Some(meth) => meth.iter().map(|m| (&m).to_string()).join(","),
            None => "NA".to_string(),
        };

        f.debug_struct("Read")
            .field("id", &self.id)
            .field("bases", &std::str::from_utf8(&self.bases).unwrap())
            .field("meth", &meth)
            .field("cigar", &self.cigar)
            .finish()
    }
}

impl HiFiRead {
    /// Creates a `HiFiRead` from an HTSlib record and a genomic region.
    ///
    /// # Arguments
    /// * `rec` - A BAM record from HTSlib.
    /// * `region` - The `GenomicRegion` struct representing the region of interest.
    ///
    /// # Returns
    /// Returns a `HiFiRead` populated with data extracted from the BAM record and the specified region.
    pub fn from_hts_rec(rec: &bam::Record, region: &GenomicRegion) -> HiFiRead {
        let id = str::from_utf8(rec.qname()).unwrap().to_string();
        let is_reverse = rec.is_reverse();
        let bases = rec.seq().as_bytes();
        let quals = rec.qual().to_vec();

        let meth = get_mm_tag(rec).and_then(|mm_tag| {
            get_ml_tag(rec)
                .and_then(|ml_tag| parse_meth_tags(mm_tag, ml_tag))
                .and_then(|tags| {
                    if rec.is_reverse() {
                        meth::decode_on_minus(&bases, &tags)
                    } else {
                        meth::decode_on_plus(&bases, &tags)
                    }
                })
        });

        let mapq = rec.mapq();
        let hp_tag = get_hp_tag(rec);
        let read_qual = get_rq_tag(rec);

        let cigar = if !rec.is_unmapped() {
            Some(Cigar {
                ref_pos: rec.reference_start(),
                ops: rec.cigar().take().to_vec(),
            })
        } else {
            None
        };

        let start_offset = (rec.reference_start() - region.start as i64) as i32;
        let end_offset = (rec.reference_end() - region.end as i64) as i32;

        let mismatch_offsets = cigar.as_ref().map(|c| snp::extract_snps_offset(c, region));

        HiFiRead {
            id,
            is_reverse,
            bases,
            quals,
            meth,
            read_qual,
            mismatch_offsets,
            start_offset,
            end_offset,
            cigar,
            hp_tag,
            mapq,
        }
    }
}

/// Parses methylation tags from a BAM record into a `MethInfo` struct.
///
/// # Arguments
/// * `mm_tag` - The MM tag from the BAM record.
/// * `ml_tag` - The ML tag from the BAM record.
///
/// # Returns
/// Returns an `Option<MethInfo>` which is `Some` if the tags could be parsed, otherwise `None`.
fn parse_meth_tags(mm_tag: Aux, ml_tag: Aux) -> Option<MethInfo> {
    let mm_tag = match mm_tag {
        Aux::String(tag) => tag,
        _ => panic!("Unexpected MM tag format: {:?}", mm_tag),
    };

    let mm_tag = mm_tag
        .strip_prefix("C+m?")
        .or_else(|| mm_tag.strip_prefix("C+m"))?;
    let mm_tag = mm_tag.trim_matches(',');

    let first_mm_rec = mm_tag.split(';').next()?;
    if first_mm_rec.is_empty() {
        return None;
    }
    let poses = first_mm_rec
        .split(',')
        .map(|n| n.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();

    let mut probs = match ml_tag {
        Aux::ArrayU8(tag) => tag.iter().collect::<Vec<_>>(),
        _ => panic!("Unexpected ML tag format: {:?}", ml_tag),
    };

    // If MM tag contains a semicolon, it must be composite
    if mm_tag.contains(';') {
        probs = probs[..poses.len()].to_vec();
    }

    assert_eq!(poses.len(), probs.len());

    Some(MethInfo { poses, probs })
}

/// Retrieves the MM tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<Aux>` which is `Some` if the MM tag is present, otherwise `None`.
fn get_mm_tag(rec: &bam::Record) -> Option<Aux> {
    rec.aux(b"MM").or_else(|_| rec.aux(b"Mm")).ok()
}

/// Retrieves the ML tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<Aux>` which is `Some` if the ML tag is present, otherwise `None`.
fn get_ml_tag(rec: &bam::Record) -> Option<Aux> {
    rec.aux(b"ML").or_else(|_| rec.aux(b"Ml")).ok()
}

/// Retrieves the RQ (read quality) tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<f64>` which is `Some` if the RQ tag is present and can be parsed as a float, otherwise `None`.
pub fn get_rq_tag(rec: &bam::Record) -> Option<f64> {
    match rec.aux(b"rq") {
        Ok(Aux::Float(value)) => Some(f64::from(value)),
        _ => None,
    }
}

/// Retrieves the HP (haplotype) tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<u8>` which is `Some` if the HP tag is present and can be parsed as a byte, otherwise `None`.
fn get_hp_tag(rec: &bam::Record) -> Option<u8> {
    match rec.aux(b"HP") {
        Ok(Aux::U8(value)) => Some(value),
        _ => None,
    }
}
