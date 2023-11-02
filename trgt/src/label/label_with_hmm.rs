use super::hmm::{self, HmmMotif};
use super::hmm_defs::HMM_DEFS;
use super::{Annotation, Span};
use crate::locus::Locus;
use itertools::Itertools;

pub fn label_with_hmm(locus: &Locus, seqs: &Vec<String>) -> Vec<Annotation> {
    let encoding = HMM_DEFS.get(&locus.struc[..]).unwrap();
    let hmm = decode_hmm(encoding);
    let mut annotations = Vec::new();

    for seq in seqs {
        let labels = hmm.label(seq);
        let labels = hmm.label_motifs(&labels);
        let motif_counts = count_motifs(&locus.motifs, &labels);
        let labels = Some(collapse_labels(labels));
        let purity = 1.0;
        annotations.push(Annotation {
            labels,
            motif_counts,
            purity,
        });
    }

    annotations
}

fn decode_hmm(encoding: &str) -> hmm::Hmm {
    let mats = encoding.split('|').collect_vec();
    assert!(mats.len() == 3);
    let ems = decode_emissions(mats[0]);
    let transitions = decode_transitions(mats[1]);

    assert!(ems.len() == transitions.len());
    let num_states = ems.len();
    let mut hmm = hmm::Hmm::new(num_states);

    for (state, state_ems) in ems.into_iter().enumerate() {
        hmm.set_ems(state, state_ems);
    }

    for (state, (in_states, probs)) in transitions.into_iter().enumerate() {
        hmm.set_trans(state, in_states, probs);
    }

    hmm.motifs = decode_motifs(mats[2]);

    hmm
}

fn decode_emissions(encoding: &str) -> Vec<Vec<f64>> {
    let mut mat = Vec::new();
    for row in encoding.split("],[") {
        let row = row.trim_matches(|c| "\"[]".contains(c));
        let row = row
            .split(',')
            .map(|e| e.parse::<f64>().unwrap())
            .collect_vec();
        mat.push(row);
    }

    mat
}

fn decode_transitions(encoding: &str) -> Vec<(Vec<usize>, Vec<f64>)> {
    let mut transitions = Vec::new();
    let encoding = encoding
        .strip_prefix("[[[")
        .unwrap()
        .strip_suffix("]]]")
        .unwrap();

    for row in encoding.split("]],[[") {
        if row.chars().all(|c| "[],".contains(c)) {
            transitions.push((Vec::new(), Vec::new()));
            continue;
        }

        let states_and_probs = row.split("],[").collect_vec();
        assert!(states_and_probs.len() == 2);
        let (states, probs) = (states_and_probs[0], states_and_probs[1]);
        let states = states
            .split(',')
            .map(|e| e.parse::<usize>().unwrap())
            .collect_vec();
        let probs = probs
            .split(',')
            .map(|e| e.parse::<f64>().unwrap())
            .collect_vec();

        assert!(states.len() == probs.len());
        transitions.push((states, probs));
    }

    transitions
}

fn decode_motifs(encoding: &str) -> Vec<HmmMotif> {
    let mut motifs = Vec::new();
    for motif_encoding in encoding.split("),(") {
        let (start_state, end_state, motif_index) = motif_encoding
            .trim_matches('(')
            .trim_end_matches(')')
            .split(',')
            .map(|m| m.parse::<usize>().unwrap())
            .collect_tuple()
            .unwrap();
        motifs.push(HmmMotif {
            start_state,
            end_state,
            motif_index,
        });
    }

    motifs
}

fn collapse_labels(spans: Vec<Span>) -> Vec<Span> {
    let mut collapsed = Vec::new();
    for span in spans {
        if collapsed.is_empty() {
            collapsed.push(span);
            continue;
        }

        let last_span = collapsed.last_mut().unwrap();
        if last_span.motif_index == span.motif_index && last_span.end == span.start {
            last_span.end = span.end;
        } else {
            collapsed.push(span);
        }
    }
    collapsed
}

fn count_motifs(motifs: &Vec<String>, labels: &Vec<Span>) -> Vec<usize> {
    let mut motif_counts = vec![0; motifs.len()];
    for span in labels {
        motif_counts[span.motif_index] += 1;
    }
    motif_counts
}
