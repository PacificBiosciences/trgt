use super::HiFiRead;
use crate::reads::cigar::{Cigar, CigarOp};

/// Clip an aligned read to a given reference region
pub fn clip_to_region(read: HiFiRead, region: (i64, i64)) -> Option<HiFiRead> {
    let cigar = read.cigar.as_ref()?;
    let (clipped_ref_start, clipped_query_start, clipped_cigar) = clip_cigar(cigar, region)?;

    let mut clipped_bases = Vec::new();
    let mut clipped_quals = Vec::new();
    let mut query_pos = clipped_query_start;

    for op in &clipped_cigar {
        let op_query_len = get_query_len(op);
        clipped_bases
            .extend(&read.bases[query_pos as usize..query_pos as usize + op_query_len as usize]);
        clipped_quals
            .extend(&read.quals[query_pos as usize..query_pos as usize + op_query_len as usize]);
        query_pos += op_query_len;
    }

    let clipped_query_end = query_pos;
    let clipped_meth = match &read.meth {
        Some(profile) => {
            let mut meth_iter = profile.iter();

            let mut clipped_meth = Vec::new();
            for index in 0..read.bases.len() - 1 {
                let dinuc = &read.bases[index..index + 2];
                let in_region =
                    clipped_query_start as usize <= index && index < clipped_query_end as usize;
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

    let cigar = Cigar {
        ref_pos: clipped_ref_start,
        ops: clipped_cigar,
    };

    Some(HiFiRead {
        bases: clipped_bases,
        quals: clipped_quals,
        meth: clipped_meth,
        cigar: Some(cigar),
        ..read
    })
}

fn get_reference_end(cigar: &Cigar) -> i64 {
    let mut ref_len = cigar.ref_pos;
    for op in &cigar.ops {
        ref_len += get_ref_len(op);
    }

    ref_len
}

/// Clips an alignment to a given reference region
///
/// Outputs reference start, query start, and operations of the clipped alignment
fn clip_cigar(cigar: &Cigar, region: (i64, i64)) -> Option<(i64, i64, Vec<CigarOp>)> {
    let (region_start, region_end) = region;
    let (read_start, read_end) = (cigar.ref_pos, get_reference_end(cigar));

    if read_end <= region_start || region_end <= read_start {
        return None;
    }

    let mut ref_pos = cigar.ref_pos;
    let mut query_pos = 0;

    let mut op_iter = cigar.ops.iter();
    let mut current_op = op_iter.next();

    let mut clipped_ops = Vec::new();

    // Skip operations outside of the target region
    while current_op.is_some() && ref_pos + get_ref_len(current_op.unwrap()) <= region_start {
        ref_pos += get_ref_len(current_op.unwrap());
        query_pos += get_query_len(current_op.unwrap());
        current_op = op_iter.next();
    }

    let mut clipped_ref_start = ref_pos;
    let mut clipped_query_start = query_pos;

    // Split operation overlapping the left flank (if any)
    if ref_pos < region_start {
        // Take care of situations where a single operation spans the entire region
        let ref_outside_len = region_start - ref_pos;
        let op_ref_len = get_ref_len(current_op.unwrap());
        let clipped_op_ref_len = if ref_pos + op_ref_len <= region_end {
            op_ref_len - ref_outside_len
        } else {
            region_end - region_start
        } as u32;
        clipped_ops.push(match current_op.unwrap() {
            CigarOp::Match(_) => CigarOp::Match(clipped_op_ref_len),
            CigarOp::RefSkip(_) => CigarOp::RefSkip(clipped_op_ref_len),
            CigarOp::Del(_) => CigarOp::Del(clipped_op_ref_len),
            CigarOp::Equal(_) => CigarOp::Equal(clipped_op_ref_len),
            CigarOp::Diff(_) => CigarOp::Diff(clipped_op_ref_len),
            op => panic!("Unexpected operation {:?}", op),
        });

        clipped_ref_start += ref_outside_len;
        if get_query_len(clipped_ops.last().unwrap()) != 0 {
            clipped_query_start += ref_outside_len;
        }

        ref_pos += get_ref_len(current_op.unwrap());
        query_pos += get_query_len(current_op.unwrap());
        current_op = op_iter.next();
    }

    // Copy operations contained within the region
    while current_op.is_some() && ref_pos + get_ref_len(current_op.unwrap()) <= region_end {
        clipped_ops.push(*current_op.unwrap());
        ref_pos += get_ref_len(current_op.unwrap());
        query_pos += get_query_len(current_op.unwrap());
        current_op = op_iter.next();
    }

    // Split operation overlapping the right flank (if any)
    if let Some(op) = current_op {
        if ref_pos < region_end {
            let ref_inside_len = (region_end - ref_pos) as u32;
            let clipped_op = match op {
                CigarOp::Match(_) => CigarOp::Match(ref_inside_len),
                CigarOp::RefSkip(_) => CigarOp::RefSkip(ref_inside_len),
                CigarOp::Del(_) => CigarOp::Del(ref_inside_len),
                CigarOp::Equal(_) => CigarOp::Equal(ref_inside_len),
                CigarOp::Diff(_) => CigarOp::Diff(ref_inside_len),
                _ => panic!("Unexpected operation {:?}", op),
            };
            clipped_ops.push(clipped_op);
        }
    }

    Some((clipped_ref_start, clipped_query_start, clipped_ops))
}

fn get_ref_len(op: &CigarOp) -> i64 {
    match op {
        CigarOp::Match(len)
        | CigarOp::RefSkip(len)
        | CigarOp::Del(len)
        | CigarOp::Equal(len)
        | CigarOp::Diff(len) => *len as i64,
        CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
    }
}

fn get_query_len(op: &CigarOp) -> i64 {
    match op {
        CigarOp::Match(len)
        | CigarOp::Equal(len)
        | CigarOp::Diff(len)
        | CigarOp::Ins(len)
        | CigarOp::SoftClip(len) => *len as i64,
        CigarOp::RefSkip(_) | CigarOp::Del(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn make_read(bases: &str, meths: Vec<u8>, cigar: Cigar) -> HiFiRead {
        HiFiRead {
            id: "read1".to_string(),
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
    fn if_no_overlap_then_none() {
        //                     /AAATC
        //           CGC--TCGTTACG
        // CCCCCCCCCCCGCGGTCATTACGCCCCCCCCCC
        // |---10---||-----13----||---10---|
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        let clipped_read = clip_to_region(read, (0, 10));
        assert_eq!(clipped_read, None);

        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        let clipped_read = clip_to_region(read, (23, 33));
        assert_eq!(clipped_read, None);
    }

    #[test]
    fn if_alignment_contained_inside_region_then_original_read() {
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        let read = make_read("AAAAACGCTCGTTAAATCACGAAAAAAAAAA", vec![10, 20, 30], cigar);
        let expected = read.clone();

        let clipped_read = clip_to_region(read, (9, 23));
        assert_eq!(clipped_read, Some(expected));
    }

    #[test]
    fn if_alignment_overlaps_left_flank_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        let clipped_read = clip_to_region(read, (0, 15));

        let cigar = make_cigar(10, "3=2D");
        let expected = make_read("CGC", vec![10], cigar);
        assert_eq!(clipped_read, Some(expected));
    }

    #[test]
    fn if_op_overlaps_flanks_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        let clipped_read = clip_to_region(read, (12, 17));

        let cigar = make_cigar(12, "1=2D2=");
        let expected = make_read("CTC", vec![20], cigar);
        assert_eq!(clipped_read, Some(expected));
    }

    #[test]
    fn if_op_spans_entire_region_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        let clipped_read = clip_to_region(read, (21, 22));

        let cigar = make_cigar(21, "1=");
        let expected = make_read("C", vec![30], cigar);
        assert_eq!(clipped_read, Some(expected));
    }

    #[test]
    fn if_alignment_starts_inside_region_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let read = make_read("CGCTCGTTAAATCACG", vec![10, 20, 30], cigar);
        let clipped_read = clip_to_region(read, (0, 17));

        let cigar = make_cigar(10, "3=2D2=");
        let expected = make_read("CGCTC", vec![10, 20], cigar);
        assert_eq!(clipped_read, Some(expected));
    }
}
