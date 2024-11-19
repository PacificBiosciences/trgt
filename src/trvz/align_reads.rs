use super::align::Align;
use super::align::{AlignOp, AlignSeg};
use super::read::{project_betas, Betas, Read};
use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment as BioAlign;
use itertools::Itertools;

type BioOp = bio::alignment::AlignmentOperation;

/// Align reads to the consensus sequence
pub fn align_reads(
    consensus: &str,
    consensus_align: &[AlignSeg],
    reads: &[&Read],
) -> Vec<(Align, Betas)> {
    let likely_max_len = consensus.len() + 50;
    let mut aligner = get_aligner(likely_max_len);
    let mut read_aligns = Vec::new();

    for read in reads {
        let bio_align = aligner.global(read.seq.as_bytes(), consensus.as_bytes());
        let align = convert(consensus_align, &bio_align);
        let betas = project_betas(&bio_align, &read.betas);
        read_aligns.push((align, betas));
    }

    read_aligns
}

/// Convert a rust-bio alignment into an internal alignment
fn convert(consensus_align: &[AlignSeg], bio_align: &BioAlign) -> Align {
    assert!(bio_align.xstart == 0);
    assert!(bio_align.ystart == 0);

    let mut seg_type_by_ref = Vec::new();
    for align_seg in consensus_align {
        let seg_type = align_seg.seg_type;
        seg_type_by_ref.extend(match align_seg.op {
            AlignOp::Del | AlignOp::Match | AlignOp::Subst => vec![seg_type; align_seg.width],
            AlignOp::Ins => vec![],
        });
    }

    let mut ref_pos = 0;
    let mut ops_and_segs = Vec::new();
    for op in &bio_align.operations {
        let seg_type = if ref_pos == seg_type_by_ref.len() {
            // Handle trailing insertion
            assert_eq!(*op, BioOp::Ins);
            seg_type_by_ref[ref_pos - 1]
        } else {
            seg_type_by_ref[ref_pos]
        };
        ops_and_segs.push((*op, seg_type));
        ref_pos += match *op {
            BioOp::Match => 1,
            BioOp::Subst => 1,
            BioOp::Del => 1,
            BioOp::Ins => 0,
            _ => panic!("Unhandled operation {op:?}"),
        };
    }

    let mut align = Vec::new();
    let mut ref_pos = 0;
    let mut seq_pos = 0;

    for ((bio_op, seg_type), group) in &ops_and_segs.iter().group_by(|rec| *rec) {
        let seg_type = *seg_type;
        let run_len = group.count();

        let align_seg = match *bio_op {
            BioOp::Match => AlignSeg {
                width: run_len,
                op: AlignOp::Match,
                seg_type,
            },
            BioOp::Subst => AlignSeg {
                width: run_len,
                op: AlignOp::Subst,
                seg_type,
            },
            BioOp::Del => AlignSeg {
                width: run_len,
                op: AlignOp::Del,
                seg_type,
            },
            BioOp::Ins => AlignSeg {
                width: 0,
                op: AlignOp::Ins,
                seg_type,
            },
            _ => panic!("No logic to handle {bio_op:?}"),
        };

        align.push(align_seg);

        ref_pos += match *bio_op {
            BioOp::Match => run_len,
            BioOp::Subst => run_len,
            BioOp::Del => run_len,
            BioOp::Ins => 0,
            _ => panic!("Unhandled operation {:?}", *bio_op),
        };

        seq_pos += match *bio_op {
            BioOp::Match => run_len,
            BioOp::Subst => run_len,
            BioOp::Del => 0,
            BioOp::Ins => run_len,
            _ => panic!("Unhandled operation {:?}", *bio_op),
        };
    }

    assert_eq!(bio_align.x_aln_len(), seq_pos);
    assert_eq!(bio_align.y_aln_len(), ref_pos);

    align
}

type ScoreFunc = fn(u8, u8) -> i32;

fn score(a: u8, b: u8) -> i32 {
    if a == b {
        1i32
    } else {
        -1i32
    }
}

fn get_aligner<'a>(likely_max_len: usize) -> Aligner<&'a ScoreFunc> {
    Aligner::with_capacity(
        likely_max_len,
        likely_max_len,
        -5,
        -1,
        &(score as ScoreFunc),
    )
}
