use super::hmm::Hmm;
use crate::cli::handle_error_and_exit;
use itertools::Itertools;

#[derive(Debug, PartialEq, Clone)]
pub enum HmmEvent {
    Match,
    Mismatch,
    Ins,
    Del,
    Trans, // Silent states that don't encode alignment operation
    Skip,  // Skip state that matches bases outside of any motif run
    MotifStart,
    MotifEnd,
}

pub fn get_events(hmm: &Hmm, motifs: &[Vec<u8>], states: &[usize], query: &[u8]) -> Vec<HmmEvent> {
    let mut state_to_hmm_motif = vec![-1; hmm.num_states];
    for hmm_motif in &hmm.motifs {
        for state in hmm_motif.start_state..=hmm_motif.end_state {
            state_to_hmm_motif[state] = hmm_motif.motif_index as i32;
        }
    }

    let mut base_index = 0;
    let mut events = Vec::new();
    let base_consumers = [
        HmmEvent::Match,
        HmmEvent::Mismatch,
        HmmEvent::Ins,
        HmmEvent::Skip,
    ];
    for state_index in 0..states.len() {
        let state = states[state_index];

        let motif_index = state_to_hmm_motif[state];
        if motif_index == -1 {
            events.push(HmmEvent::Trans);
            continue;
        }

        let hmm_motif = &hmm.motifs[motif_index as usize];

        if state == hmm_motif.start_state {
            events.push(HmmEvent::MotifStart);
            // This is fine since the last state must have motif_index of -1
            let next_state = states[state_index + 1];
            let num_dels = next_state - state - 1;
            events.extend(vec![HmmEvent::Del; num_dels]);
            continue;
        }

        if state == hmm_motif.end_state {
            events.push(HmmEvent::MotifEnd);
            continue;
        }

        if motif_index as usize + 1 == hmm.motifs.len() {
            events.push(HmmEvent::Skip);
            base_index += 1;
            continue;
        }

        let offset = state - hmm_motif.start_state - 1;
        let motif_len = motifs[hmm_motif.motif_index].len();

        let event = match offset.div_euclid(motif_len) {
            0 => {
                let base = query[base_index];
                if base == get_base_match(&hmm, state) {
                    HmmEvent::Match
                } else {
                    HmmEvent::Mismatch
                }
            }
            1 => HmmEvent::Ins,
            2 => HmmEvent::Del,
            _ => handle_error_and_exit(format!("Event decoding error")),
        };
        if base_consumers.contains(&event) {
            base_index += 1;
        }
        events.push(event);
    }

    events
}

pub fn get_base_match(hmm: &Hmm, state: usize) -> u8 {
    let ems = &hmm.ems[state];
    assert_eq!(ems.len(), 5);

    if !hmm.emits_base(state) {
        return b' ';
    }

    let max_lp = *ems.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let top_indexes = ems
        .iter()
        .enumerate()
        .filter(|(_i, p)| **p == max_lp)
        .map(|(i, _p)| i)
        .collect_vec();
    if top_indexes.len() == 1 {
        match top_indexes[0] {
            0 => b'#',
            1 => b'A',
            2 => b'T',
            3 => b'C',
            4 => b'G',
            _ => handle_error_and_exit(format!("Unexpected base match event")),
        }
    } else {
        b' '
    }
}

#[cfg(test)]
mod tests {
    use super::super::builder::build_hmm;
    use super::*;

    #[test]
    fn states_match_the_most_likely_base() {
        let motifs = vec!["A".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        assert_eq!(get_base_match(&hmm, 3), b'A');
    }

    #[test]
    fn silent_states_match_a_blank_character() {
        let mut hmm = Hmm::new(3);
        hmm.set_ems(0, vec![1.00, 0.00, 0.00, 0.00, 0.00]); // Start state
        hmm.set_ems(1, vec![0.00, 0.50, 0.50, 0.00, 0.00]); // State with identical scores for 'A' and 'T'
        hmm.set_ems(2, vec![1.00, 0.00, 0.00, 0.00, 0.00]); // End state
        assert_eq!(get_base_match(&hmm, 1), b' ',);
    }
}
