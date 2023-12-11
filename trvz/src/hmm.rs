use itertools::Itertools;
use std::collections::HashMap;

// List of abbreviations
// lp = log probability
// ems = emissions

#[derive(Debug, Clone, PartialEq)]
pub struct Span {
    pub motif_index: usize,
    pub start: usize,
    pub end: usize,
}

type MatF64 = Vec<Vec<f64>>;
type MatInt = Vec<Vec<usize>>;

#[derive(Debug, PartialEq)]
pub struct Hmm {
    num_states: usize,
    ems: MatF64,
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
        let mut ems = Vec::with_capacity(num_states);
        for _ in 0..num_states {
            ems.push(vec![f64::NEG_INFINITY; 5]);
        }

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

    fn calc_viterbi_score(
        &self,
        query: &[u8],
        scores: &MatF64,
        state: usize,
        index: usize,
    ) -> Option<(usize, f64)> {
        let symbol = query[index];
        let is_silent = self.ems[state].iter().all(|e| e.is_infinite());
        let em_term = if is_silent {
            0.0
        } else {
            self.ems[state][symbol as usize]
        };
        let lookback = if is_silent { 0 } else { 1 };

        let in_states = &self.in_states[state];

        let mut max_score = f64::NEG_INFINITY;
        let mut best_state = None;

        if index == 0 && !in_states.is_empty() && lookback == 1 {
            return None;
        }

        for (in_state_index, prev_state) in in_states.iter().enumerate() {
            let prev_score = scores[*prev_state][index - lookback];
            let trans_lp = self.in_lps[state][in_state_index];
            let in_state_score = prev_score + trans_lp + em_term;

            if in_state_score > max_score {
                best_state = Some(prev_state);
                max_score = in_state_score
            }
        }

        //TODO: Investigate "em_term != 0.0"
        if index == 0 && in_states.is_empty() && em_term.is_finite() {
            max_score = em_term;
            best_state = Some(&state);
        }

        best_state.map(|state| (*state, max_score))
    }

    fn generate_mats(&self, query: &Vec<u8>) -> (Vec<Vec<f64>>, Vec<Vec<Option<usize>>>) {
        let ordered_states = self.order_states();
        let mut scores = vec![vec![f64::NEG_INFINITY; query.len()]; self.num_states];
        let mut states: Vec<Vec<Option<usize>>> = vec![vec![None; query.len()]; self.num_states];

        for index in 0..query.len() {
            for state in &ordered_states {
                let result = self.calc_viterbi_score(query, &scores, *state, index);
                if let Some((prev_state, score)) = result {
                    scores[*state][index] = score;
                    states[*state][index] = Some(prev_state);
                }
            }
        }
        (scores, states)
    }

    /*
    fn summarize(&self, traceback_states: &Vec<usize>) -> Vec<(usize, usize)> {
        let mut count_by_state = Vec::new();
        for (state, group) in &traceback_states.iter().group_by(|v| **v) {
            count_by_state.push((state, group.count()));
        }
        count_by_state
    } */

    fn traceback(&self, query: &Vec<u8>, states: &[Vec<Option<usize>>]) -> Vec<usize> {
        let mut state = self.num_states - 1;
        let exp_max_capacity = self.num_states + (self.num_states as f64 * 0.20).round() as usize;
        let mut traceback_states = Vec::with_capacity(exp_max_capacity);
        let mut index = query.len() - 1;
        while state != 0 {
            traceback_states.push(state);

            let prev_state = states[state][index].unwrap();
            if self.ems[state].iter().any(|e| e.is_finite()) {
                index -= 1;
            }
            state = prev_state;
        }
        traceback_states.push(0);
        traceback_states.reverse();
        traceback_states
    }

    pub fn label(&self, query: &str) -> Vec<usize> {
        let query = "#"
            .bytes()
            .chain(query.bytes().chain("#".bytes()))
            .map(encode_base)
            .collect_vec();

        let (_scores, states) = self.generate_mats(&query);
        self.traceback(&query, &states)
    }

    pub fn label_motifs(&self, states: &Vec<usize>) -> Vec<Span> {
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

            if state_to_motif.contains_key(&state) {
                let motif_index = state_to_motif[&state];
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

                let motif_start = if motif_spans.is_empty() {
                    0
                } else {
                    motif_spans.last().unwrap().end
                };
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

    fn emits_base(&self, state: usize) -> bool {
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
