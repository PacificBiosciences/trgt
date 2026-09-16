use super::spans::Span;
use itertools::Itertools;
use std::collections::HashMap;

// List of abbreviations
// lp = log probability
// ems = emissions

type MatF64 = Vec<Vec<f64>>;
type MatInt = Vec<Vec<usize>>;

#[derive(Debug, PartialEq)]
pub struct Hmm {
    pub num_states: usize,
    pub ems: MatF64,
    in_states: MatInt,
    in_lps: MatF64,
    pub motifs: Vec<HmmMotif>,
}

#[derive(Debug, PartialEq)]
pub struct HmmMotif {
    pub start_state: usize,
    pub end_state: usize,
    pub motif_index: usize,
}

impl Hmm {
    pub fn new(num_states: usize) -> Hmm {
        let ems = vec![vec![f64::NEG_INFINITY; 5]; num_states];
        let in_states = vec![Vec::new(); num_states];
        let in_lps = vec![Vec::new(); num_states];
        let motifs = Vec::new();

        Hmm {
            num_states,
            ems,
            in_states,
            in_lps,
            motifs,
        }
    }

    pub fn set_trans(&mut self, target_state: usize, in_states: Vec<usize>, in_probs: Vec<f64>) {
        self.in_states[target_state] = in_states;
        self.in_lps[target_state] = in_probs.iter().map(|v| v.ln()).collect_vec();
    }

    pub fn set_ems(&mut self, target_state: usize, ems: Vec<f64>) {
        assert!(ems.len() == 5 || ems.is_empty());
        self.ems[target_state] = ems.iter().map(|v| v.ln()).collect_vec();
    }

    fn generate_traceback(&self, query: &[u8]) -> Vec<usize> {
        let ordered_states = self.order_states();
        let silent: Vec<bool> = self
            .ems
            .iter()
            .map(|ems| ems.iter().all(|e| e.is_infinite()))
            .collect();
        // Emitting states read the previous column; silent states read the current
        // column in topological order. Earlier score columns are never needed.
        let mut previous = vec![f64::NEG_INFINITY; self.num_states];
        let mut current = vec![f64::NEG_INFINITY; self.num_states];
        // Store predecessors by query position for contiguous writes. MAX marks
        // unreachable cells without the extra word used by Option<usize>.
        let mut states = vec![usize::MAX; query.len() * self.num_states];

        for (index, &symbol) in query.iter().enumerate() {
            current.fill(f64::NEG_INFINITY);
            let column = &mut states[index * self.num_states..(index + 1) * self.num_states];
            for &state in &ordered_states {
                let in_states = &self.in_states[state];
                if index == 0 && !in_states.is_empty() && !silent[state] {
                    continue;
                }
                let em_term = if silent[state] {
                    0.0
                } else {
                    self.ems[state][symbol as usize]
                };
                let scores = if silent[state] { &current } else { &previous };
                let mut max_score = f64::NEG_INFINITY;
                let mut best_state = usize::MAX;
                for (in_state_index, &prev_state) in in_states.iter().enumerate() {
                    // Preserve arithmetic and predecessor order, including ties.
                    let score = scores[prev_state] + self.in_lps[state][in_state_index] + em_term;
                    if score > max_score {
                        best_state = prev_state;
                        max_score = score;
                    }
                }
                if index == 0 && in_states.is_empty() && em_term.is_finite() {
                    max_score = em_term;
                    best_state = state;
                }
                current[state] = max_score;
                column[state] = best_state;
            }
            std::mem::swap(&mut previous, &mut current);
        }
        states
    }

    /*
    fn summarize(&self, traceback_states: &Vec<usize>) -> Vec<(usize, usize)> {
        let mut count_by_state = Vec::new();
        for (state, group) in &traceback_states.iter().group_by(|v| **v) {
            count_by_state.push((state, group.count()));
        }
        count_by_state
    } */

    fn traceback(&self, query: &[u8], states: &[usize]) -> Vec<usize> {
        let mut state = self.num_states - 1;
        let exp_max_capacity = self.num_states + (self.num_states as f64 * 0.20).round() as usize;
        let mut traceback_states = Vec::with_capacity(exp_max_capacity);
        let mut index = query.len() - 1;
        while state != 0 {
            traceback_states.push(state);

            let prev_state = states[index * self.num_states + state];
            assert_ne!(prev_state, usize::MAX, "Unreachable traceback state");
            if self.ems[state].iter().any(|e| e.is_finite()) {
                index -= 1;
            }
            state = prev_state;
        }
        traceback_states.push(0);
        traceback_states.reverse();
        traceback_states
    }

    pub fn label(&self, query: &[u8]) -> Vec<usize> {
        if query.is_empty() {
            return Vec::new();
        }
        let query = b"#"
            .iter()
            .chain(query.iter().chain(b"#".iter()))
            .map(|&b| encode_base(b))
            .collect_vec();

        let states = self.generate_traceback(&query);
        self.traceback(&query, &states)
    }

    pub fn label_motifs(&self, states: &[usize]) -> Vec<Span> {
        let state_to_motif: HashMap<usize, usize> = self
            .motifs
            .iter()
            .enumerate()
            .map(|(index, m)| (m.start_state, index))
            .collect();

        let mut motif_spans: Vec<Span> = Vec::new();
        let mut state_index = 0;

        while state_index < states.len() {
            let state = states[state_index];
            if let Some(&motif_index) = state_to_motif.get(&state) {
                let motif = &self.motifs[motif_index];
                let mut motif_span = 0;

                while states[state_index] != motif.end_state {
                    motif_span += self.emits_base(states[state_index]) as usize;
                    state_index += 1;
                }

                while states[state_index] == motif.end_state {
                    motif_span += self.emits_base(states[state_index]) as usize;
                    state_index += 1;
                }

                let motif_start = motif_spans.last().map_or(0, |span| span.end);
                let motif_end = motif_start + motif_span;

                motif_spans.push(Span {
                    motif_index,
                    start: motif_start,
                    end: motif_end,
                });
            } else {
                assert!(!self.emits_base(state));
                state_index += 1;
            }
        }

        motif_spans
    }

    pub fn emits_base(&self, state: usize) -> bool {
        self.ems[state].iter().skip(1).any(|e| e.is_finite())
    }

    fn order_states(&self) -> Vec<usize> {
        let mut normal_states = Vec::new();
        let mut silent_states = Vec::new();

        for state in 0..self.num_states {
            if self.ems[state].iter().all(|e| e.is_infinite()) {
                silent_states.push(state);
            } else {
                normal_states.push(state);
            }
        }

        let mut sorted = Vec::new();

        while !silent_states.is_empty() {
            let mut unused = Vec::new();

            for state in &silent_states {
                let has_incoming_silent = self.in_states[*state]
                    .iter()
                    .any(|s| silent_states.contains(s));
                if !has_incoming_silent {
                    sorted.push(*state)
                } else {
                    unused.push(*state);
                }
            }

            assert!(unused.len() < silent_states.len());
            silent_states = unused;
        }

        normal_states.extend(sorted);
        normal_states
    }
}

fn encode_base(base: u8) -> u8 {
    match base {
        b'#' => 0,
        b'A' => 1,
        b'T' => 2,
        b'C' => 3,
        b'G' => 4,
        _ => panic!(),
    }
}

#[cfg(test)]
mod equivalence_tests {
    use super::*;
    use crate::hmm::build_hmm;
    use rand::{Rng, SeedableRng};

    // Independent full-matrix recurrence from the original implementation.
    fn reference_label(hmm: &Hmm, query: &[u8]) -> Vec<usize> {
        if query.is_empty() {
            return Vec::new();
        }
        let query: Vec<u8> = b"#"
            .iter()
            .chain(query.iter())
            .chain(b"#".iter())
            .map(|&b| encode_base(b))
            .collect();
        let mut scores = vec![vec![f64::NEG_INFINITY; query.len()]; hmm.num_states];
        let mut states = vec![vec![None; query.len()]; hmm.num_states];
        let ordered_states = hmm.order_states();
        for index in 0..query.len() {
            for &state in &ordered_states {
                let silent = hmm.ems[state].iter().all(|e| e.is_infinite());
                let em_term = if silent {
                    0.0
                } else {
                    hmm.ems[state][query[index] as usize]
                };
                let lookback = !silent as usize;
                let incoming = &hmm.in_states[state];
                if index == 0 && !incoming.is_empty() && lookback == 1 {
                    continue;
                }
                let mut max_score = f64::NEG_INFINITY;
                let mut best = None;
                for (i, &prev) in incoming.iter().enumerate() {
                    let score = scores[prev][index - lookback] + hmm.in_lps[state][i] + em_term;
                    if score > max_score {
                        max_score = score;
                        best = Some(prev);
                    }
                }
                if index == 0 && incoming.is_empty() && em_term.is_finite() {
                    max_score = em_term;
                    best = Some(state);
                }
                if best.is_some() {
                    scores[state][index] = max_score;
                    states[state][index] = best;
                }
            }
        }
        let mut path = Vec::new();
        let mut state = hmm.num_states - 1;
        let mut index = query.len() - 1;
        while state != 0 {
            path.push(state);
            let prev = states[state][index].unwrap();
            if hmm.ems[state].iter().any(|e| e.is_finite()) {
                index -= 1;
            }
            state = prev;
        }
        path.push(0);
        path.reverse();
        path
    }

    #[test]
    fn rolling_scores_match_full_matrix_exhaustively() {
        for motifs in [
            vec![b"A".to_vec()],
            vec![b"CAG".to_vec()],
            vec![b"AC".to_vec(), b"A".to_vec()],
            vec![b"N".to_vec(), b"AN".to_vec()],
            vec![b"CAG".to_vec(), b"CAG".to_vec()], // tied paths
            vec![b"ACGTACGT".to_vec()],             // silent deletion chains
        ] {
            let hmm = build_hmm(&motifs);
            for len in 0..=5 {
                for mut encoded in 0..4_usize.pow(len) {
                    let query: Vec<u8> = (0..len)
                        .map(|_| {
                            let base = b"ACGT"[encoded % 4];
                            encoded /= 4;
                            base
                        })
                        .collect();
                    assert_eq!(
                        hmm.label(&query),
                        reference_label(&hmm, &query),
                        "{motifs:?} {query:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn rolling_scores_match_full_matrix_for_long_and_mutated_repeats() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1729);
        for case in 0..200 {
            let motifs: Vec<Vec<u8>> = (0..rng.random_range(1..=4))
                .map(|_| {
                    (0..rng.random_range(1..=40))
                        .map(|_| b"ACGTN"[rng.random_range(0..5)])
                        .collect()
                })
                .collect();
            let mut query = Vec::new();
            for _ in 0..rng.random_range(1..=20) {
                for &base in &motifs[rng.random_range(0..motifs.len())] {
                    if rng.random_range(0..10) == 0 {
                        continue;
                    }
                    query.push(if base == b'N' || rng.random_range(0..10) == 0 {
                        b"ACGT"[rng.random_range(0..4)]
                    } else {
                        base
                    });
                    if rng.random_range(0..10) == 0 {
                        query.push(b"ACGT"[rng.random_range(0..4)]);
                    }
                }
            }
            let hmm = build_hmm(&motifs);
            assert_eq!(
                hmm.label(&query),
                reference_label(&hmm, &query),
                "case {case}"
            );
        }
        for len in [1, 6, 100] {
            let motif: Vec<u8> = (0..len).map(|i| b"ACGT"[i % 4]).collect();
            let hmm = build_hmm(&[motif.clone()]);
            let query = motif.repeat(3);
            assert_eq!(hmm.label(&query), reference_label(&hmm, &query));
        }
    }

    // A 1000-base motif builds ~3000 states, where reference_label's full
    // num_states x query_len matrices dominate the whole suite's runtime. Real
    // catalog motifs top out at 24 bases, so this stress case runs on demand:
    // cargo test -- --ignored
    #[test]
    #[ignore]
    fn rolling_scores_match_full_matrix_for_very_long_motif() {
        let motif: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
        let hmm = build_hmm(&[motif.clone()]);
        let query = motif.repeat(3);
        assert_eq!(hmm.label(&query), reference_label(&hmm, &query));
    }
}

/*#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn make_ct_hmm() -> Hmm {
        let mut hmm = Hmm::new(5);
        //                                    #     A     T     C     G
        hmm.set_ems(0, vec![1.00, 0.00, 0.00, 0.00, 0.00]);
        hmm.set_ems(1, vec![0.00, 0.25, 0.25, 0.25, 0.25]);
        hmm.set_ems(2, vec![0.00, 0.10, 0.40, 0.40, 0.10]);
        hmm.set_ems(3, vec![0.00, 0.25, 0.25, 0.25, 0.25]);
        hmm.set_ems(4, vec![1.00, 0.00, 0.00, 0.00, 0.00]);

        hmm.set_trans(1, vec![0, 1], vec![0.5, 0.5]);
        hmm.set_trans(2, vec![0, 1, 2], vec![0.5, 0.5, 0.4]);
        hmm.set_trans(3, vec![2, 3], vec![0.3, 0.5]);
        hmm.set_trans(4, vec![2, 3], vec![0.3, 0.5]);

        hmm
    }

    #[test]
    fn create_hmm_by_setting_state_number() {
        let hmm = Hmm::new(5);
        assert_eq!(hmm.ems.len(), hmm.num_states);
    }

    #[test]
    fn set_emissions_with_corresponding_fn() {
        let mut hmm = Hmm::new(5);
        let state = 0;
        //                                        #    A    T    C    G
        hmm.set_ems(state, vec![1.0, 0.0, 0.0, 0.0, 0.0]);
        assert_relative_eq!(hmm.ems[state][0], 1.0_f64.ln());
        assert_relative_eq!(hmm.ems[state][1], 0.0_f64.ln());
        assert_relative_eq!(hmm.ems[state][2], 0.0_f64.ln());
        assert_relative_eq!(hmm.ems[state][3], 0.0_f64.ln());
        assert_relative_eq!(hmm.ems[state][4], 0.0_f64.ln());
    }

    #[test]
    fn set_transitions_with_corresponding_fn() {
        let mut hmm = Hmm::new(5);
        let target_state = 2;
        let in_states = vec![0, 1];
        let in_probs = vec![0.5, 0.5];
        hmm.set_trans(target_state, in_states, in_probs);

        assert_eq!(hmm.in_states[target_state][0], 0);
        assert_eq!(hmm.in_states[target_state][1], 1);
        assert_relative_eq!(hmm.in_lps[target_state][0], 0.5_f64.ln());
        assert_relative_eq!(hmm.in_lps[target_state][1], 0.5_f64.ln());
    }

    #[test]
    fn use_hmm_to_find_ct_runs() {
        let hmm = make_ct_hmm();
        let query = "AAAATCTCTCTCGGGG";
        let labels = hmm.label(query);

        let expected = vec![(0, 1), (1, 4), (2, 8), (3, 4), (4, 1)];
        assert_eq!(hmm.summarize(&labels), expected);
    }

    #[test]
    fn use_silent_states_to_model_deletions() {
        let mut hmm = Hmm::new(8);
        //                                    #     A     T     C     G
        hmm.set_ems(0, vec![1.00, 0.00, 0.00, 0.00, 0.00]);
        hmm.set_ems(1, vec![0.00, 0.70, 0.10, 0.10, 0.10]);
        hmm.set_ems(2, vec![0.00, 0.10, 0.10, 0.70, 0.10]);
        hmm.set_ems(3, vec![0.00, 0.10, 0.70, 0.10, 0.10]);
        hmm.set_ems(4, vec![0.00, 0.70, 0.10, 0.10, 0.10]);
        hmm.set_ems(5, vec![]);
        hmm.set_ems(6, vec![]);
        hmm.set_ems(7, vec![1.00, 0.00, 0.00, 0.00, 0.00]);

        hmm.set_trans(1, vec![0, 1], vec![1.0, 0.4]);
        hmm.set_trans(2, vec![1], vec![0.3]);
        hmm.set_trans(3, vec![2, 5], vec![0.5, 0.5]);
        hmm.set_trans(4, vec![3, 4, 6], vec![1.0, 0.5, 1.0]);
        hmm.set_trans(5, vec![1], vec![0.3]);
        hmm.set_trans(6, vec![5, 2], vec![0.5, 0.5]);
        hmm.set_trans(7, vec![4], vec![0.5]);

        let expected = vec![(0, 1), (1, 1), (5, 1), (6, 1), (4, 7), (7, 1)];
        assert_eq!(hmm.summarize(&hmm.label("AAAAAAAA")), expected);

        let expected = vec![(0, 1), (1, 4), (2, 1), (3, 1), (4, 3), (7, 1)];
        assert_eq!(hmm.summarize(&hmm.label("AAAACTAAA")), expected);

        let expected = vec![(0, 1), (1, 2), (5, 1), (3, 1), (4, 1), (7, 1)];
        assert_eq!(hmm.summarize(&hmm.label("AATA")), expected);
    }

    #[test]
    fn use_hmms_with_loops_to_find_motif_copies() {
        //  0   5    1   2   3   4    6
        //  B---RS---A---T---T---RE---E
        //      |________________|
        //

        let mut hmm = Hmm::new(7);
        //                                    #     A     T     C     G
        hmm.set_ems(0, vec![1.00, 0.00, 0.00, 0.00, 0.00]);
        hmm.set_ems(1, vec![0.00, 0.70, 0.10, 0.10, 0.10]);
        hmm.set_ems(2, vec![0.00, 0.10, 0.70, 0.10, 0.10]);
        hmm.set_ems(3, vec![0.00, 0.10, 0.70, 0.10, 0.10]);
        hmm.set_ems(4, vec![]);
        hmm.set_ems(5, vec![]);
        hmm.set_ems(6, vec![1.00, 0.00, 0.00, 0.00, 0.00]);

        hmm.set_trans(1, vec![5], vec![1.0]);
        hmm.set_trans(2, vec![1], vec![1.0]);
        hmm.set_trans(3, vec![2], vec![1.0]);
        hmm.set_trans(4, vec![3], vec![1.0]);
        hmm.set_trans(5, vec![0, 4], vec![1.0, 0.5]);
        hmm.set_trans(6, vec![4], vec![0.5]);

        let labels = hmm.label("ATTATTATTATTATT");
        let motifs = hmm.label_motifs(&labels, (4, 5));
        assert_eq!(motifs, vec![(0, 3), (3, 6), (6, 9), (9, 12), (12, 15)]);
    }
}
*/
