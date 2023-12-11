use super::struc::*;
use bio::alignment::distance::simd::*;

pub fn refine_motif_counts(
    motifs: &[Motif],
    query: &str,
    motif_counts: &[usize],
) -> (Vec<usize>, f64) {
    let mut top_counts = apply_multiplicity_counts(motifs, motif_counts);
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

fn get_candidates(motifs: &[Motif], seed_counts: &[usize]) -> Vec<Vec<usize>> {
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

fn get_top_candidate(motifs: &[Motif], query: &str, mut candidates: Vec<Vec<usize>>) -> Vec<usize> {
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

fn apply_multiplicity_counts(motifs: &[Motif], motif_counts: &[usize]) -> Vec<usize> {
    motifs
        .iter()
        .zip(motif_counts.iter())
        .map(
            |(motif, &count)| {
                if motif.mult == Mult::Many {
                    count
                } else {
                    1
                }
            },
        )
        .collect()
}

fn get_dist(template: &[Motif], query: &str, motif_counts: &[usize]) -> u32 {
    let mut reference = String::new();
    for (motif, &count) in template.iter().zip(motif_counts) {
        reference.push_str(&motif.seq.repeat(count));
    }
    levenshtein(reference.as_bytes(), query.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_motifs() -> Vec<Motif> {
        vec![
            Motif {
                seq: "A".to_string(),
                mult: Mult::Once,
            },
            Motif {
                seq: "C".to_string(),
                mult: Mult::Many,
            },
            Motif {
                seq: "G".to_string(),
                mult: Mult::Once,
            },
        ]
    }

    #[test]
    fn test_get_dist_exact_match() {
        let template = setup_motifs();
        let query = "ACG";
        let motif_counts = vec![1, 1, 1];
        assert_eq!(get_dist(&template, query, &motif_counts), 0);
    }

    #[test]
    fn test_get_dist_1_ins() {
        let template = setup_motifs();
        let query = "AACG";
        let motif_counts = vec![1, 1, 1];
        assert_eq!(get_dist(&template, query, &motif_counts), 1);
    }

    #[test]
    fn test_get_dist_1_del() {
        let template = setup_motifs();
        let query = "CG";
        let motif_counts = vec![1, 1, 1];
        assert_eq!(get_dist(&template, query, &motif_counts), 1);
    }

    #[test]
    fn test_get_dist_1_sub() {
        let template = setup_motifs();
        let query = "ACC";
        let motif_counts = vec![1, 1, 1];
        assert_eq!(get_dist(&template, query, &motif_counts), 1);
    }

    #[test]
    fn test_get_dist_diff() {
        let template = setup_motifs();
        let query = "TTT";
        let motif_counts = vec![1, 1, 1];
        assert_eq!(get_dist(&template, query, &motif_counts), 3);
    }

    #[test]
    fn test_get_dist_empty_template() {
        let template: Vec<Motif> = Vec::new();
        let query = "ACG";
        let motif_counts: Vec<usize> = Vec::new();
        assert_eq!(
            get_dist(&template, query, &motif_counts),
            query.len() as u32
        );
    }

    #[test]
    fn test_get_dist_empty_query() {
        let template = setup_motifs();
        let query = "";
        let motif_counts = vec![1, 2, 1];
        // Distance should be the sum of motif_counts since query is empty.
        assert_eq!(
            get_dist(&template, query, &motif_counts),
            motif_counts.iter().sum::<usize>() as u32
        );
    }

    #[test]
    fn test_get_dist_empty_both() {
        let template: Vec<Motif> = Vec::new();
        let query = "";
        let motif_counts: Vec<usize> = Vec::new();
        assert_eq!(get_dist(&template, query, &motif_counts), 0);
    }

    #[test]
    fn test_apply_multiplicity_counts_basic() {
        let motifs = vec![
            Motif {
                seq: "A".to_string(),
                mult: Mult::Once,
            },
            Motif {
                seq: "B".to_string(),
                mult: Mult::Many,
            },
        ];
        let motif_counts = vec![1, 3];
        let expected = vec![1, 3];
        assert_eq!(apply_multiplicity_counts(&motifs, &motif_counts), expected);
    }

    #[test]
    fn test_apply_multiplicity_counts_empty_vectors() {
        let motifs: Vec<Motif> = Vec::new();
        let motif_counts: Vec<usize> = Vec::new();
        let expected: Vec<usize> = Vec::new();
        assert_eq!(apply_multiplicity_counts(&motifs, &motif_counts), expected);
    }

    #[test]
    fn test_all_once() {
        let motifs = vec![
            Motif {
                seq: "A".to_string(),
                mult: Mult::Once,
            },
            Motif {
                seq: "B".to_string(),
                mult: Mult::Once,
            },
        ];
        let motif_counts = vec![5, 10]; // These counts should be ignored.
        let expected = vec![1, 1];
        assert_eq!(apply_multiplicity_counts(&motifs, &motif_counts), expected);
    }

    #[test]
    fn test_all_many() {
        let motifs = vec![
            Motif {
                seq: "A".to_string(),
                mult: Mult::Many,
            },
            Motif {
                seq: "B".to_string(),
                mult: Mult::Many,
            },
        ];
        let motif_counts = vec![2, 4];
        let expected = vec![2, 4];
        assert_eq!(apply_multiplicity_counts(&motifs, &motif_counts), expected);
    }

    #[test]
    fn test_mixed_multiplicities() {
        let motifs = vec![
            Motif {
                seq: "A".to_string(),
                mult: Mult::Once,
            },
            Motif {
                seq: "B".to_string(),
                mult: Mult::Many,
            },
            Motif {
                seq: "C".to_string(),
                mult: Mult::Once,
            },
            Motif {
                seq: "D".to_string(),
                mult: Mult::Many,
            },
        ];
        let motif_counts = vec![5, 10, 15, 20];
        let expected = vec![1, 10, 1, 20];
        assert_eq!(apply_multiplicity_counts(&motifs, &motif_counts), expected);
    }
}
