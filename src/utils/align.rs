use crate::{
    commands::genotype::THREAD_WFA_CONSENSUS,
    wfaligner::{CigarOp, WFAligner},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Clone, Copy)]
pub struct TrgtScoring {
    pub mism_scr: i32,
    pub gapo_scr: i32,
    pub gape_scr: i32,
}

pub fn align(backbone: &str, seqs: &[&str]) -> Vec<Vec<CigarOp>> {
    let backbone_bytes = backbone.as_bytes();
    let alignments: Vec<Vec<CigarOp>> = seqs
        .par_iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            THREAD_WFA_CONSENSUS.with(|aligner_cell| {
                let mut aligner = aligner_cell.borrow_mut();
                let _status = aligner.align_end_to_end(backbone_bytes, seq_bytes);
                WFAligner::decode_sam_cigar(&aligner.get_sam_cigar(true))
            })
        })
        .collect();
    alignments
}
