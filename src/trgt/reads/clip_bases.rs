use super::{
    cigar::{Cigar, CigarOp, CigarOpExt},
    HiFiRead,
};
use crate::utils::handle_error_and_exit;

/// Clip a specified number of bases from both sides of the readx
impl HiFiRead {
    pub fn clip_bases(&self, left_len: usize, right_len: usize) -> Option<HiFiRead> {
        if left_len + right_len >= self.bases.len() {
            return None;
        }

        let clipped_bases = self.bases[left_len..self.bases.len() - right_len].to_vec();
        let clipped_quals = self.quals[left_len..self.quals.len() - right_len].to_vec();

        let clipped_cigar = if self.cigar.is_some() {
            clip_cigar(self.cigar.as_ref().unwrap().clone(), left_len, right_len)
        } else {
            None
        };

        let clipped_meth = match &self.meth {
            Some(profile) => {
                let (region_start, region_end) = (left_len, self.bases.len() - right_len);
                let mut meth_iter = profile.iter();

                let mut clipped_meth = Vec::new();
                for index in 0..self.bases.len() - 1 {
                    let dinuc = &self.bases[index..index + 2];
                    let in_region = region_start <= index && index < region_end;
                    if dinuc == b"CG" {
                        if in_region {
                            clipped_meth.push(*meth_iter.next().unwrap());
                        } else {
                            meth_iter.next();
                        }
                    }
                }

                Some(clipped_meth)
            }
            None => None,
        };

        Some(HiFiRead {
            bases: clipped_bases,
            quals: clipped_quals,
            meth: clipped_meth,
            cigar: clipped_cigar,
            id: self.id.clone(),
            mismatch_offsets: self.mismatch_offsets.clone(),
            ..*self
        })
    }
}

/// Clips CIGAR by a specified number of bases on left and right
fn clip_cigar(cigar: Cigar, mut left_len: usize, right_len: usize) -> Option<Cigar> {
    let align_query_len: usize = cigar.ops.iter().map(|op| op.get_query_len() as usize).sum();
    assert!(align_query_len >= left_len + right_len);
    let mut keep_len = (align_query_len - left_len - right_len) as u32;

    let mut op_iter = cigar.ops.into_iter();
    let mut current_op = op_iter.next();
    let mut ref_pos: usize = cigar.ref_pos as usize;

    while left_len != 0 {
        let query_len = current_op.unwrap().get_query_len() as usize;
        if query_len > left_len {
            // Split the operation
            let leftover_len = (query_len - left_len) as u32;
            current_op = match current_op.unwrap() {
                CigarOp::Match(_) => Some(CigarOp::Match(leftover_len)),
                CigarOp::Diff(_) => Some(CigarOp::Diff(leftover_len)),
                CigarOp::Ins(_) => Some(CigarOp::Ins(leftover_len)),
                CigarOp::Equal(_) => Some(CigarOp::Equal(leftover_len)),
                CigarOp::SoftClip(_) => Some(CigarOp::SoftClip(leftover_len)),
                // Since op has non-zero query length, it can't be Del, RefSkip, HardClip, or Pad
                _ => handle_error_and_exit(format!("Unexpected operation {current_op:?}")),
            };

            if current_op.unwrap().get_ref_len() != 0 {
                ref_pos += left_len
            };
            left_len = 0;
        } else {
            left_len -= query_len;
            ref_pos += current_op.unwrap().get_ref_len() as usize;
            current_op = op_iter.next();
        }
    }

    let mut clipped_ops = Vec::new();
    while current_op.is_some() && keep_len != 0 {
        let query_len = current_op.unwrap().get_query_len() as u32;
        if query_len > keep_len {
            clipped_ops.push(match current_op.unwrap() {
                CigarOp::Match(_) => CigarOp::Match(keep_len),
                CigarOp::Diff(_) => CigarOp::Diff(keep_len),
                CigarOp::Ins(_) => CigarOp::Ins(keep_len),
                CigarOp::Equal(_) => CigarOp::Equal(keep_len),
                CigarOp::SoftClip(_) => CigarOp::SoftClip(keep_len),
                // Since op has non-zero query length, it can't be Del, RefSkip, HardClip, or Pad
                _ => handle_error_and_exit(format!("Unexpected operation {current_op:?}")),
            });
            keep_len = 0;
        } else {
            keep_len -= query_len;
            clipped_ops.push(current_op.unwrap());
            current_op = op_iter.next();
        }
    }

    Some(Cigar {
        ref_pos: ref_pos as i64,
        ops: clipped_ops,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn make_read(bases: &str, meths: Vec<u8>, cigar: Cigar) -> HiFiRead {
        HiFiRead {
            id: "read".to_string(),
            is_reverse: false,
            bases: bases.as_bytes().to_vec(),
            quals: "(".repeat(bases.len()).as_bytes().to_vec(),
            meth: Some(meths),
            read_qual: None,
            mismatch_offsets: None,
            start_offset: 0,
            end_offset: 0,
            cigar: Some(cigar),
            hp_tag: None,
            mapq: 60,
        }
    }

    fn make_cigar(ref_pos: i64, encoding: &str) -> Cigar {
        let ops = CigarString::try_from(encoding).unwrap().to_vec();
        Cigar { ref_pos, ops }
    }

    #[test]
    fn get_none_if_clip_whole_query() {
        //                     /AAATC
        //           CGC--TCGTTACG
        // CCCCCCCCCCCGCGGTCATTACGCCCCCCCCCC
        // |---10---||-----13----||---10---|
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);

        assert_eq!(read.clip_bases(16, 0), None);
        assert_eq!(read.clip_bases(0, 16), None);
        assert_eq!(read.clip_bases(12, 4), None);
    }

    #[test]
    fn clip_from_left() {
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        let read = make_read("AAAAACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);

        let cigar = make_cigar(10, "2S3=2D2=1X2=5I3=10S");
        let expected = make_read("AACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);
        assert_eq!(read.clip_bases(3, 0), Some(expected));

        let cigar = make_cigar(10, "3=2D2=1X2=5I3=10S");
        let expected = make_read("CGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);
        assert_eq!(read.clip_bases(5, 0), Some(expected));

        let cigar = make_cigar(17, "1X2=5I3=10S");
        let expected = make_read("GTTAAATCACGAAAAAAAAAA", vec![30], cigar);
        assert_eq!(read.clip_bases(10, 0), Some(expected));
    }

    #[test]
    fn clip_from_right() {
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        let read = make_read("AAAAACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);

        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=5S");
        let expected = make_read("AAAAACGCTCGTTAAATCACGAAAAA", vec![10, 20, 30], cigar);
        assert_eq!(read.clip_bases(0, 5), Some(expected));

        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=");
        let expected = make_read("AAAAACGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        assert_eq!(read.clip_bases(0, 10), Some(expected));

        let cigar = make_cigar(10, "5S3=2D2=1X2=3I");
        let expected = make_read("AAAAACGCTCGTTAAA", vec![10, 20], cigar);
        assert_eq!(read.clip_bases(0, 15), Some(expected));
    }

    #[test]
    fn clip_from_both_sides() {
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        let read = make_read("AAAAACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);

        let cigar = make_cigar(10, "3=2D2=1X2=5I3=5S");
        let expected = make_read("CGCTCGTTAAATCACGAAAAA", vec![10, 20, 30], cigar);
        assert_eq!(read.clip_bases(5, 5), Some(expected));

        let cigar = make_cigar(13, "2D2=1X2=5I2=");
        let expected = make_read("TCGTTAAATCAC", vec![20, 30], cigar);
        assert_eq!(read.clip_bases(8, 11), Some(expected));
    }

    #[test]
    fn get_none_if_clip_all_bases() {
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        let read = make_read("AAAAACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);

        assert_eq!(read.clip_bases(31, 0), None);
        assert_eq!(read.clip_bases(0, 31), None);
        assert_eq!(read.clip_bases(30, 30), None);
    }

    #[test]
    fn cigar_generated_even_if_clipped_alignment_has_no_query_len() {
        // It might be better to mark the read as unaligned in such cases (especially if IGV complains)
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        let read = make_read("AAAAACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);

        let cigar = make_cigar(20, "5I");
        let expected = make_read("AAATC", Vec::new(), cigar);
        assert_eq!(read.clip_bases(13, 13), Some(expected));
    }
}
