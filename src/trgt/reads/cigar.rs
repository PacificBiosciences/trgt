pub type CigarOp = rust_htslib::bam::record::Cigar;

pub trait CigarOpExt {
    fn get_ref_len(&self) -> i64;
    fn get_query_len(&self) -> i64;
}

impl CigarOpExt for CigarOp {
    fn get_ref_len(&self) -> i64 {
        match self {
            CigarOp::Match(len)
            | CigarOp::RefSkip(len)
            | CigarOp::Del(len)
            | CigarOp::Equal(len)
            | CigarOp::Diff(len) => *len as i64,
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
        }
    }

    fn get_query_len(&self) -> i64 {
        match self {
            CigarOp::Match(len)
            | CigarOp::Equal(len)
            | CigarOp::Diff(len)
            | CigarOp::Ins(len)
            | CigarOp::SoftClip(len) => *len as i64,
            CigarOp::RefSkip(_) | CigarOp::Del(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Cigar {
    pub ref_pos: i64,
    pub ops: Vec<CigarOp>,
}

impl Cigar {
    #[allow(dead_code)]
    pub fn query_len(&self) -> usize {
        self.ops.iter().map(|op| op.get_query_len() as usize).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_len() {
        let cigar = Cigar {
            ref_pos: 0,
            ops: vec![
                CigarOp::Match(10),
                CigarOp::Ins(5),
                CigarOp::Del(3),
                CigarOp::SoftClip(2),
            ],
        };
        assert_eq!(cigar.query_len(), 17);
    }

    #[test]
    fn test_get_ref_len() {
        assert_eq!(CigarOp::Match(10).get_ref_len(), 10);
        assert_eq!(CigarOp::Ins(5).get_ref_len(), 0);
        assert_eq!(CigarOp::Del(3).get_ref_len(), 3);
        assert_eq!(CigarOp::SoftClip(2).get_ref_len(), 0);
    }

    #[test]
    fn test_get_query_len() {
        assert_eq!(CigarOp::Match(10).get_query_len(), 10);
        assert_eq!(CigarOp::Ins(5).get_query_len(), 5);
        assert_eq!(CigarOp::Del(3).get_query_len(), 0);
        assert_eq!(CigarOp::SoftClip(2).get_query_len(), 2);
    }
}
