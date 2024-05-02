use bio::alignment::{pairwise::Aligner, Alignment};
use itertools::Itertools;

#[derive(Debug, Clone, Copy)]
pub struct TrgtScoring {
    pub match_scr: i32,
    pub mism_scr: i32,
    pub gapo_scr: i32,
    pub gape_scr: i32,
    pub kmer_len: usize,
    pub bandwidth: usize,
}

pub fn align(backbone: &str, seqs: &[&str]) -> Vec<Alignment> {
    let mut aligner = Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });
    seqs.iter()
        .map(|seq| aligner.global(seq.as_bytes(), backbone.as_bytes()))
        .collect_vec()
}
