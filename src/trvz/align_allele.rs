use super::align::AlleleAlign;
use super::align_consensus::align_consensus;
use super::align_reads::align_reads;
use super::locus::Locus;
use super::read::Read;

pub fn get_allele_align(locus: &Locus, consensus: &str, reads: &[&Read]) -> AlleleAlign {
    let (consensus_align, motif_bounds) = align_consensus(locus, consensus);
    let read_aligns = align_reads(consensus, &consensus_align, reads);
    AlleleAlign {
        seq: (consensus_align, motif_bounds),
        reads: read_aligns,
    }
}
