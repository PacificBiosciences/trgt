//! Module for representing and building read information from alignment records.
//!

use super::{cigar::Cigar, snp};
use crate::utils::GenomicRegion;
use itertools::Itertools;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::str;

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

fn get_meth(rec: &bam::Record, bases: &[u8]) -> Option<Vec<u8>> {
    let reverse = rec.is_reverse();
    let cpg_indices: Vec<usize> = std::str::from_utf8(bases)
        .unwrap()
        .match_indices("CG")
        .map(|(x, _)| {
            x + if reverse { 1 } else { 0 } // need Gs for reverse
        })
        .collect::<Vec<usize>>();
    let num_cpgs = cpg_indices.len();
    let mut ans: Vec<u8> = vec![0; cpg_indices.len()];
    let mut ind = 0;
    if let Ok(mods) = rec.basemods_iter() {
        for (pos, m) in mods.flatten() {
            if m.canonical_base as u8 == b'C' {
                let pos = pos as usize;
                while ind < num_cpgs && cpg_indices[ind] < pos {
                    ind += 1;
                }
                if ind < num_cpgs && pos == cpg_indices[ind] {
                    ans[ind] = m.qual as u8;
                    ind += 1;
                }
            }
        }
        if ind == 0 {
            //empty MM/ML
            return None;
        }
        if reverse {
            ans.reverse();
        }
        return Some(ans);
    }
    None
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

        let meth: Option<Vec<u8>> = get_meth(rec, &bases);

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
#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{record::Aux, Record};

    fn create_record(bases: &[u8], mm: &str, ml: &[u8], reverse: bool) -> Record {
        let mut rec = Record::new();
        let qual = vec![40; bases.len()];
        if reverse {
            rec.set_flags(0x10);
        }
        rec.set(b"test", None, bases, &qual);
        rec.push_aux(b"MM", Aux::String(mm)).unwrap();
        rec.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
        rec
    }

    #[test]
    fn test_basemods_error() {
        let bases = b"ACGTCG";
        let rec = create_record(bases, "no", &[42], false);
        let res = get_meth(&rec, bases);
        assert!(res.is_none());
    }

    #[test]
    fn test_matching_modifications() {
        let bases = b"AGTCTAGACTCCGTAATTACTCGCCTAG";
        let mm = "C+m,3,1;";
        let ml = [249, 4];
        let rec = create_record(bases, mm, &ml, false);
        let res = get_meth(&rec, bases);
        assert_eq!(res, Some(vec![249, 4]));
    }
}
