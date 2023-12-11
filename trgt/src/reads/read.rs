use super::{cigar::Cigar, meth, snp::extract_snps_offset};
use crate::utils::GenomicRegion;
use itertools::Itertools;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux};
use std::str;

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
    pub hp_tag: Option<u8>,
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
    pub fn from_hts_rec(rec: bam::Record, region: &GenomicRegion) -> HiFiRead {
        let id = str::from_utf8(rec.qname()).unwrap().to_string();
        let bases = rec.seq().as_bytes();

        let meth = get_mm_tag(&rec).and_then(|mm_tag| {
            get_ml_tag(&rec)
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
        let hp_tag = get_hp_tag(&rec);
        let read_qual = get_rq_tag(&rec);

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

        let mismatch_offsets = cigar.as_ref().map(|c| extract_snps_offset(c, region));

        HiFiRead {
            id,
            bases,
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
    rec.aux(b"MM").or_else(|_| rec.aux(b"Mm")).ok()
}

fn get_ml_tag(rec: &bam::Record) -> Option<Aux> {
    rec.aux(b"ML").or_else(|_| rec.aux(b"Ml")).ok()
}

fn get_rq_tag(rec: &bam::Record) -> Option<f64> {
    match rec.aux(b"rq") {
        Ok(Aux::Float(value)) => Some(f64::from(value)),
        _ => None,
    }
}

fn get_hp_tag(rec: &bam::Record) -> Option<u8> {
    match rec.aux(b"HP") {
        Ok(Aux::U8(value)) => Some(u8::from(value)),
        _ => None,
    }
}
