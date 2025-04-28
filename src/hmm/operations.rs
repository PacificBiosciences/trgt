use super::Hmm;
use crate::hmm::hmm_model::HmmMotif;
use std::collections::{HashMap, HashSet};

/// Replace imperfect STR motif copies with Skip states
pub fn remove_imperfect_motifs(
    hmm: &Hmm,
    motifs: &[Vec<u8>],
    states: &[usize],
    query: &[u8],
    max_motif_len: usize,
) -> Vec<usize> {
    if states.is_empty() {
        return Vec::new();
    }
    let start_state_to_motif: HashMap<usize, &HmmMotif> =
        hmm.motifs.iter().map(|m| (m.start_state, m)).collect();

    assert!(states.len() > 4);
    let mut updated_states = vec![states[0], states[1]];

    let motif_start_states: HashSet<usize> = hmm.motifs.iter().map(|m| m.start_state).collect();
    let motif_end_states: HashSet<usize> = hmm.motifs.iter().map(|m| m.end_state).collect();
    let motif_run_end_state = hmm.num_states - 2;

    let mut state_index = 2;
    let mut base_index = 0;

    while state_index != states.len() {
        assert!(motif_start_states.contains(&states[state_index]));
        let mut motif_states = Vec::new();
        let mut motif_sequence = String::new();
        while !motif_end_states.contains(&states[state_index]) {
            motif_states.push(states[state_index]);
            if hmm.emits_base(states[state_index]) {
                motif_sequence.push(query[base_index] as char); // get_base_match(hmm, states[state_index])
                base_index += 1;
            }
            state_index += 1;
        }
        motif_states.push(states[state_index]);
        state_index += 1;

        let motif_rec = start_state_to_motif[motif_states.first().unwrap()];
        let motif_len = (motif_rec.end_state - motif_rec.start_state) / 3;
        let mut keep_motif_match = true;
        // Only STR motif matches can be removed
        let skip_motif = motif_rec.motif_index + 1 == hmm.motifs.len();
        if !skip_motif && motif_len <= max_motif_len {
            let motif = &motifs[motif_rec.motif_index];
            if motif_sequence.len() < motif.len() {
                keep_motif_match = false;
            } else {
                let iterator = motif.iter().zip(motif_sequence.chars().take(motif.len()));
                for (expected, observed) in iterator {
                    if *expected as char != 'N' && (observed != *expected as char) {
                        keep_motif_match = false;
                    }
                }
            }
        }

        if keep_motif_match {
            updated_states.extend(motif_states.iter());
        } else {
            let bases_consumed = motif_states.iter().filter(|s| hmm.emits_base(**s)).count();
            let skip_motif = hmm.motifs.last().unwrap();
            updated_states.push(skip_motif.start_state);
            updated_states.extend(vec![skip_motif.start_state + 1; bases_consumed]);
            updated_states.push(skip_motif.end_state);
        }

        if states[state_index] == motif_run_end_state {
            updated_states.extend(&states[state_index..state_index + 2]);
            state_index += 2;
        }
    }

    updated_states
}
