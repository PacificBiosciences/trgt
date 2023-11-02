use super::{struc::decode_regexp, Annotation, Span};
use bio::alignment::distance::simd::*;

pub fn label_motif(struc: &str, seqs: &Vec<String>) -> Vec<Annotation> {
    let motif = decode_regexp(struc).pop().unwrap().seq;
    let mut annotations = Vec::new();

    for seq in seqs {
        let labels = Some(vec![Span {
            motif_index: 0,
            start: 0,
            end: seq.len(),
        }]);
        let purity = calc_purity(&motif, seq);
        let motif_counts = vec![(seq.len() as f64 / motif.len() as f64).round() as usize];
        annotations.push(Annotation {
            labels,
            motif_counts,
            purity,
        });
    }

    annotations
}

fn calc_purity(motif: &str, seq: &str) -> f64 {
    // Assume that alleles of zero length are 100% pure
    if seq.is_empty() {
        return 1.0;
    }
    let num_motifs = seq.len() / motif.len();
    let leftover_len = seq.len() % motif.len();
    let pure_tr = motif.repeat(num_motifs) + &motif[..leftover_len];
    let dist = levenshtein(pure_tr.as_bytes(), seq.as_bytes());
    let similarity = seq.len() - dist as usize;
    similarity as f64 / seq.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_calc_purity_pure() {
        assert_eq!(calc_purity(&"ATGC", &"ATGCATGCATGCATGCATG"), 1.0);
    }
}
