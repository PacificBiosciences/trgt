use super::hmm::{self, HmmMotif};
use super::hmm_defs::HMM_DEFS;
use crate::locus::{BaseLabel, Locus};
use itertools::Itertools;

pub fn label_with_hmm(locus: &Locus, alleles: &Vec<String>) -> Vec<Vec<BaseLabel>> {
    let mut labels_by_allele = Vec::new();
    let encoding = HMM_DEFS.get(&locus.struc[..]).unwrap();
    let hmm = decode_hmm(encoding);
    for allele in alleles {
        let query = &allele[locus.left_flank.len()..allele.len() - locus.right_flank.len()];
        let mut labels = vec![BaseLabel::Match; locus.left_flank.len()];
        let states = hmm.label(query);
        labels.extend(get_base_labels(&hmm, &states));
        labels.extend(vec![BaseLabel::Match; locus.right_flank.len()]);
        labels_by_allele.push(labels);
    }

    labels_by_allele
}

fn get_base_labels(hmm: &hmm::Hmm, states: &Vec<usize>) -> Vec<BaseLabel> {
    let mut base_labels = vec![BaseLabel::MotifBound];
    for (start, end, _motif_index) in hmm.label_motifs(states) {
        base_labels.extend(vec![BaseLabel::Match; end - start]);
        base_labels.push(BaseLabel::MotifBound);
    }
    base_labels
}

/*pub fn label_hmm(locus: &Locus, consensuses: &Vec<String>) -> Vec<Option<Spans>> {
    let hmm = decode_hmm(&locus.struc, &locus.motifs);
    let mut labels_by_hap = Vec::new();

    for seq in consensuses {
        let labels = hmm.label(seq);
        let spans = hmm.label_motifs(&labels);
        labels_by_hap.push(Some(spans));
    }

    labels_by_hap
} */

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
