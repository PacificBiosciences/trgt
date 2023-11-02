use super::struc::*;
use bio::alignment::distance::simd::*;

pub fn refine_motif_counts(
    motifs: &Vec<Motif>,
    query: &str,
    motif_counts: &[usize],
) -> (Vec<usize>, f64) {
    let mut top_counts = fix(motifs, motif_counts);
    let mut min_dist = get_dist(motifs, query, &top_counts);

    loop {
        let candidates = get_candidates(motifs, &top_counts);
        let top_candidate = get_top_candidate(motifs, query, candidates);
        let candidate_dist = get_dist(motifs, query, &top_candidate);

        if candidate_dist < min_dist {
            min_dist = candidate_dist;
            top_counts = top_candidate;
        } else {
            break;
        }
    }

    let purity = (query.len() as f64 - min_dist as f64) / query.len() as f64;
    (top_counts, purity)
}

fn get_candidates(motifs: &Vec<Motif>, seed_counts: &[usize]) -> Vec<Vec<usize>> {
    let mut candidates = Vec::new();
    for motif_index in 0..motifs.len() {
        if motifs[motif_index].mult == Mult::Once {
            continue;
        }
        candidates.push(seed_counts.to_owned());
        candidates.last_mut().unwrap()[motif_index] += 1;

        if seed_counts[motif_index] != 0 {
            candidates.push(seed_counts.to_owned());
            candidates.last_mut().unwrap()[motif_index] -= 1;
        }
    }

    candidates
}

fn get_top_candidate(
    motifs: &Vec<Motif>,
    query: &str,
    mut candidates: Vec<Vec<usize>>,
) -> Vec<usize> {
    let mut top_candidate = candidates.pop().unwrap();
    let mut top_dist = get_dist(motifs, query, &top_candidate);
    for candidate in candidates {
        let candidate_dist = get_dist(motifs, query, &candidate);
        if candidate_dist <= top_dist {
            top_candidate = candidate;
            top_dist = candidate_dist;
        }
    }
    top_candidate
}

fn fix(motifs: &Vec<Motif>, motif_counts: &[usize]) -> Vec<usize> {
    let mut fixed_counts = Vec::new();

    for index in 0..motifs.len() {
        if motifs[index].mult == Mult::Many {
            fixed_counts.push(motif_counts[index]);
        } else {
            fixed_counts.push(1);
        }
    }

    fixed_counts
}

fn get_dist(template: &Vec<Motif>, query: &str, motif_counts: &[usize]) -> u32 {
    let mut reference = String::new();
    for index in 0..template.len() {
        let seq = &template[index];
        reference += &seq.seq.repeat(motif_counts[index]);
    }

    let score = levenshtein(reference.as_bytes(), query.as_bytes());
    score
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use crate::label::struc::decode_regexp;

    #[test]
    fn test_imperfect_htt_repeat() {
        let motifs = decode_regexp("(CAG)nCAACAG(CCG)n");
        //                 111111111111111333333
        let query = "CAGCAGCAGCATCAGCCGCCG";
        let counts = vec![5, 0, 2];
        let expected = vec![
            Span {
                motif_index: 0,
                start: 0,
                end: 9,
            },
            Span {
                motif_index: 1,
                start: 9,
                end: 15,
            },
            Span {
                motif_index: 2,
                start: 15,
                end: 21,
            },
        ];
        assert_eq!(refine_motif_counts(&motifs, &query, &counts), expected);

        //                 111111111111111333333
        //let query = "CACCAGCAGCATCAGCGGCCG";
        //let counts = vec![1, 0, 5];
        //assert_eq!(match_with_align(&motifs, &query, &counts), vec![3, 2, 2]);
    }
}

*/
