use super::struc::Motif;
pub struct KmerFilter {
    motifs: Vec<String>,
}

impl KmerFilter {
    pub fn new(motifs: &Vec<Motif>) -> KmerFilter {
        let mut kmers = Vec::new();
        for motif in motifs {
            kmers.push(motif.seq.clone());
        }
        KmerFilter { motifs: kmers }
    }

    pub fn is_passing(&self, query: &str) -> bool {
        let mut match_span = 0;
        for index in 0..query.len() {
            for motif in &self.motifs {
                if index + motif.len() >= query.len() {
                    continue;
                }

                let kmer = &query[index..index + motif.len()];
                if motif == kmer {
                    match_span += motif.len();
                }
            }
        }

        let span_frac = match_span as f64 / query.len() as f64;
        span_frac > 0.5
    }
}

#[cfg(test)]
mod tests {
    use crate::label::kmer_filter::KmerFilter;
    use crate::label::struc::decode_regexp;

    #[test]
    fn test_htt_repeat() {
        let motifs = decode_regexp("(CAG)nCCACCG(CCG)n");
        let filter = KmerFilter::new(&motifs);
        assert!(filter.is_passing("CCGCCGCCGCCG"));
        assert!(!filter.is_passing("CGGCGGCGGCGG"));
    }
}
