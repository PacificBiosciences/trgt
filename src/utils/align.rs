use crate::{commands::genotype::THREAD_WFA_CONSENSUS, wfa_aligner::cigar_wfa_to_ops};
use bio::alignment::AlignmentOperation;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Clone, Copy)]
pub struct TrgtScoring {
    pub match_scr: i32,
    pub mism_scr: i32,
    pub gapo_scr: i32,
    pub gape_scr: i32,
}

pub fn align(backbone: &str, seqs: &[&str]) -> Vec<Vec<AlignmentOperation>> {
    let backbone_bytes = backbone.as_bytes();
    let alignments: Vec<Vec<AlignmentOperation>> = seqs
        .par_iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            THREAD_WFA_CONSENSUS.with(|aligner_cell| {
                let mut aligner = aligner_cell.borrow_mut();
                let _status = aligner.align_end_to_end(backbone_bytes, seq_bytes);
                cigar_wfa_to_ops(&aligner.cigar_wfa())
            })
        })
        .collect();
    alignments
}
