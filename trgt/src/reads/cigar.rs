pub type CigarOp = rust_htslib::bam::record::Cigar;

#[derive(Debug, PartialEq, Clone)]
pub struct Cigar {
    pub ref_pos: i64,
    pub ops: Vec<CigarOp>,
}

impl Cigar {
    #[allow(dead_code)]
    pub fn query_len(&self) -> usize {
        self.ops.iter().map(|op| get_query_len(op) as usize).sum()
    }
}

pub fn get_ref_len(op: &CigarOp) -> i64 {
    match op {
        CigarOp::Match(len)
        | CigarOp::RefSkip(len)
        | CigarOp::Del(len)
        | CigarOp::Equal(len)
        | CigarOp::Diff(len) => *len as i64,
        CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
    }
}

pub fn get_query_len(op: &CigarOp) -> i64 {
    match op {
        CigarOp::Match(len)
        | CigarOp::Equal(len)
        | CigarOp::Diff(len)
        | CigarOp::Ins(len)
        | CigarOp::SoftClip(len) => *len as i64,
        CigarOp::RefSkip(_) | CigarOp::Del(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
    }
}
