use super::{hmm_model::HmmMotif, Hmm};
use itertools::Itertools;

pub fn build_hmm(motifs: &[Vec<u8>]) -> Hmm {
    // 2 terminal states + 2 run start states + 3 states of the skip block + (4n - 1) states for each motif of length n
    let num_states = 7 + motifs.iter().map(|m| 3 * m.len() + 1).sum::<usize>();
    let mut hmm = Hmm::new(num_states);

    let start = 0;
    let end = num_states - 1;
    let rs = start + 1;
    let re = end - 1;

    //                                        #     A     T     C     G
    hmm.set_ems(start, vec![1.00, 0.00, 0.00, 0.00, 0.00]);
    hmm.set_ems(end, vec![1.00, 0.00, 0.00, 0.00, 0.00]);
    hmm.set_trans(end, vec![re], vec![0.10]);

    hmm.set_ems(rs, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
    hmm.set_trans(rs, vec![start, re], vec![1.00, 0.90]);

    let rs_to_ms = 0.50; // / (motifs.len() as f64 + 1.0); <- No longer an HMM because of this change
    let me_to_re = 0.05;
    let mut mes = Vec::with_capacity(motifs.len() + 1);
    let mut ms = rs + 1;
    for motif in motifs {
        let num_motif_states = 3 * motif.len() + 1;
        let me = ms + num_motif_states - 1;

        hmm.set_ems(ms, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
        hmm.set_trans(ms, vec![rs, me], vec![rs_to_ms, 1.0 - me_to_re]);

        define_motif_block(&mut hmm, ms, motif);

        mes.push(me);
        ms += num_motif_states;
    }

    assert_eq!(ms + 3, re);

    // Defined the skip block
    let (skip_state, me) = (ms + 1, ms + 2);
    hmm.set_ems(ms, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
    hmm.set_trans(ms, vec![rs, me], vec![rs_to_ms, 1.0 - me_to_re]);

    let skip_to_skip = 0.9;
    hmm.set_ems(skip_state, vec![0.00, 0.25, 0.25, 0.25, 0.25]);
    hmm.set_trans(skip_state, vec![ms, skip_state], vec![1.0, skip_to_skip]);

    hmm.set_ems(me, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
    hmm.set_trans(me, vec![skip_state], vec![1.0 - skip_to_skip]);

    mes.push(me);

    // Define the re state
    hmm.set_ems(re, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
    hmm.set_trans(re, mes.clone(), vec![me_to_re; motifs.len() + 1]);

    // Define motif spans
    for (motif_index, motif) in motifs.iter().enumerate() {
        let me = mes[motif_index];
        let ms = me - 3 * motif.len();
        hmm.motifs.push(HmmMotif {
            start_state: ms,
            end_state: me,
            motif_index,
        });
    }

    // Add skip state span
    hmm.motifs.push(HmmMotif {
        start_state: skip_state - 1,
        end_state: skip_state + 1,
        motif_index: motifs.len(),
    });

    hmm
}

fn define_motif_block(hmm: &mut Hmm, ms: usize, motif: &[u8]) {
    let match_states = (ms + 1..ms + 1 + motif.len()).collect_vec();
    let first_ins_state = *match_states.last().unwrap() + 1;
    let ins_states = (first_ins_state..first_ins_state + motif.len()).collect_vec();
    let first_del_state = *ins_states.last().unwrap() + 1; // If any
    let del_states = (first_del_state..first_del_state + motif.len() - 1).collect_vec();

    let match_prob = 0.90;
    let ins_to_ins = 0.25;
    let match_to_indel = (1.00 - match_prob) / 2.00;
    let del_to_match = 0.50;

    // Define match states
    let mismatch_seed_prob = 2.00 * (1.00 - match_prob) / (motif.len() * (motif.len() - 1)) as f64;
    for (match_index, match_state) in match_states.iter().enumerate() {
        hmm.set_ems(*match_state, get_match_emissions(motif[match_index]));
        if match_index == 0 {
            hmm.set_trans(*match_state, vec![ms], vec![match_prob]);
        } else if match_index == 1 {
            let multiplier = motif.len() - match_index;
            let mismatch_prob = mismatch_seed_prob * multiplier as f64;
            let prev_ins = ins_states[match_index - 1];

            hmm.set_trans(
                *match_state,
                vec![match_state - 1, ms, prev_ins],
                vec![match_prob, mismatch_prob, 1.0 - ins_to_ins],
            );
        } else {
            let multiplier = motif.len() - match_index;
            let mismatch_prob = mismatch_seed_prob * multiplier as f64;
            let prev_ins = ins_states[match_index - 1];
            let prev_del = del_states[match_index - 2];

            hmm.set_trans(
                *match_state,
                vec![match_state - 1, ms, prev_ins, prev_del],
                vec![match_prob, mismatch_prob, 1.0 - ins_to_ins, del_to_match],
            );
        }
    }

    // Define insertion states
    for (ins_index, ins_state) in ins_states.iter().enumerate() {
        hmm.set_ems(*ins_state, vec![0.00, 0.25, 0.25, 0.25, 0.25]);
        let match_state = match_states[ins_index];
        hmm.set_trans(
            *ins_state,
            vec![*ins_state, match_state],
            vec![ins_to_ins, match_to_indel],
        );
    }

    // Define deletion states
    for (del_index, del_state) in del_states.iter().enumerate() {
        hmm.set_ems(*del_state, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
        let prev_match = match_states[del_index];
        if del_index == 0 {
            hmm.set_trans(*del_state, vec![prev_match], vec![match_to_indel]);
        } else {
            let prev_del = del_states[del_index - 1];
            hmm.set_trans(
                *del_state,
                vec![prev_match, prev_del],
                vec![match_to_indel, 1.0 - del_to_match],
            );
        }
    }

    let num_motif_states = 3 * motif.len() + 1;
    let me = ms + num_motif_states - 1;
    hmm.set_ems(me, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
    if !del_states.is_empty() {
        let last_match = *match_states.last().unwrap();
        let last_ins = *ins_states.last().unwrap();
        let last_del = *del_states.last().unwrap();
        hmm.set_trans(
            me,
            vec![last_match, last_ins, last_del],
            vec![match_prob, 1.0 - ins_to_ins, 1.0],
        );
    } else if !ins_states.is_empty() {
        let last_match = *match_states.last().unwrap();
        let last_ins = *ins_states.last().unwrap();
        hmm.set_trans(
            me,
            vec![last_match, last_ins],
            vec![match_prob, 1.0 - ins_to_ins],
        );
    } else {
        let last_match = *match_states.last().unwrap();
        hmm.set_trans(me, vec![last_match], vec![match_prob]);
    }
}

fn get_match_emissions(base: u8) -> Vec<f64> {
    match base {
        b'A' => vec![0.00, 0.90, 0.03, 0.03, 0.03],
        b'T' => vec![0.00, 0.03, 0.90, 0.03, 0.03],
        b'C' => vec![0.00, 0.03, 0.03, 0.90, 0.03],
        b'G' => vec![0.00, 0.03, 0.03, 0.03, 0.90],
        b'N' => vec![0.00, 0.25, 0.25, 0.25, 0.25],
        _ => panic!("Encountered unknown base {}", base as char),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmm::spans::Span;

    fn summarize(spans: &[Span]) -> Vec<(usize, usize, usize)> {
        let mut summary = Vec::new();
        for (motif_index, group) in &spans
            .iter()
            .map(|s| (s.start, s.end, s.motif_index))
            .group_by(|(_s, _e, m)| *m)
        {
            let group = group.collect_vec();
            summary.push((
                group.first().unwrap().0,
                group.last().unwrap().1,
                motif_index,
            ));
        }
        summary
    }

    #[test]
    fn annotate_two_perfect_motif_runs() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "A".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let labels = hmm.label_motifs(&hmm.label("CAGCAGCAGCAGAAAAA"));
        let expected = vec![(0, 12, 0), (12, 17, 1)];

        assert_eq!(summarize(&labels), expected);
    }

    #[test]
    fn annotate_motif_runs_separated_by_insertion() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "A".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let labels = hmm.label_motifs(&hmm.label("CAGCAGATCGATCGATCGATCGAAAAA"));
        let expected = vec![(0, 6, 0), (6, 22, 2), (22, 27, 1)];

        assert_eq!(summarize(&labels), expected);
    }

    #[test]
    fn annotate_imperfect_repeat_run() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "A".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let labels = hmm.label_motifs(&hmm.label("CAGCAGCTGCAGCAGAAACAG"));
        let expected = vec![(0, 21, 0)];

        assert_eq!(summarize(&labels), expected);
    }

    #[test]
    fn parse_aga_repeat() {
        // TODO: Consider improving this segmentation
        let motifs = vec!["AAG".as_bytes().to_vec(), "CAAC".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCCAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAATAGAAATGTGTTTAAGAATTCCTCAATAAG";
        let labels = hmm.label_motifs(&hmm.label(query));
        let expected = vec![
            (0, 6, 2),
            (6, 14, 1),
            (14, 36, 2),
            (36, 101, 0),
            (101, 119, 2),
            (119, 125, 0),
        ];

        assert_eq!(summarize(&labels), expected);
    }
}
