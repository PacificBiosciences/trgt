use super::cigar::Cigar;
use crate::utils::GenomicRegion;
use itertools::Itertools;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::str;

use super::meth;

#[derive(Debug)]
pub struct MethInfo {
    pub poses: Vec<usize>,
    pub probs: Vec<u8>,
}

#[derive(PartialEq, Clone)]
pub struct HiFiRead {
    pub id: String,
    pub bases: Vec<u8>,
    pub meth: Option<Vec<u8>>,
    pub read_qual: Option<f64>,
    pub mismatch_offsets: Option<Vec<i32>>,
    pub start_offset: i32,
    pub end_offset: i32,
    pub cigar: Option<Cigar>,
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
    pub fn from_hts_rec(rec: bam::Record, region: &GenomicRegion) -> HiFiRead {
        let id = str::from_utf8(rec.qname()).unwrap().to_string();
        let bases = rec.seq().as_bytes();

        let mm_tag = get_mm_tag(&rec);
        let ml_tag = get_ml_tag(&rec);

        let meth = match mm_tag {
            Some(mm_tag) => parse_meth_tags(mm_tag, ml_tag.unwrap()),
            None => None,
        };

        let meth = match meth {
            Some(tags) => {
                if rec.is_reverse() {
                    meth::decode_on_minus(&bases, &tags)
                } else {
                    meth::decode_on_plus(&bases, &tags)
                }
            }
            None => None,
        };

        let read_qual = get_rq_tag(&rec);

        let cigar = if rec.is_unmapped() {
            None
        } else {
            let ref_pos = rec.reference_start();
            let ops = rec.cigar().take().to_vec();
            let cigar = Cigar { ref_pos, ops };
            Some(cigar)
        };

        let start_offset = (rec.reference_start() - region.start as i64) as i32;
        let end_offset = (rec.reference_end() - region.end as i64) as i32;

        HiFiRead {
            id,
            bases,
            meth,
            read_qual,
            mismatch_offsets: None,
            start_offset,
            end_offset,
            cigar,
        }
    }
}

fn parse_meth_tags(mm_tag: Aux, ml_tag: Aux) -> Option<MethInfo> {
    let mm_tag = match mm_tag {
        Aux::String(tag) => tag,
        _ => panic!("Unexpected MM tag format: {:?}", mm_tag),
    };

    if mm_tag.is_empty() || mm_tag.len() <= 5 {
        return None;
    }

    let mm_tag = if mm_tag.contains('?') {
        &mm_tag[5..mm_tag.len() - 1]
    } else {
        &mm_tag[4..mm_tag.len() - 1]
    };
    let mm_tag = mm_tag
        .split(',')
        .map(|n| n.parse::<usize>().unwrap())
        .collect::<Vec<_>>();

    let ml_tag = match ml_tag {
        Aux::ArrayU8(tag) => tag.iter().collect::<Vec<_>>(),
        _ => panic!("Unexpected ML tag format: {:?}", ml_tag),
    };

    assert_eq!(mm_tag.len(), ml_tag.len());

    Some(MethInfo {
        poses: mm_tag,
        probs: ml_tag,
    })
}

fn get_mm_tag(rec: &bam::Record) -> Option<Aux> {
    if let Ok(value) = rec.aux(b"MM") {
        Some(value)
    } else if let Ok(value) = rec.aux(b"Mm") {
        Some(value)
    } else {
        None
    }
}

fn get_ml_tag(rec: &bam::Record) -> Option<Aux> {
    if let Ok(value) = rec.aux(b"ML") {
        Some(value)
    } else if let Ok(value) = rec.aux(b"Ml") {
        Some(value)
    } else {
        None
    }
}

fn get_rq_tag(rec: &bam::Record) -> Option<f64> {
    let rq_tag = rec.aux(b"rq");
    if rq_tag.is_err() {
        return None;
    }

    let rq_tag = rq_tag.unwrap();
    if let Aux::Float(value) = rq_tag {
        return Some(value as f64);
    }

    panic!("Unexpected rq tag format: {:?}", rq_tag);
}
