use super::hmm::{Hmm, HmmMotif};
use super::spans::Span;
use super::Annotation;
use crate::cli::handle_error_and_exit;
use crate::locus::Locus;
use itertools::Itertools;

pub fn label_with_hmm(locus: &Locus, seqs: &Vec<String>) -> Vec<Annotation> {
    let motifs = locus
        .motifs
        .iter()
        .map(|m| m.as_bytes().to_vec())
        .collect_vec();
    let hmm = build_hmm(&motifs);

    let mut annotations = Vec::new();
    for seq in seqs {
        let seq: String = seq
            .as_bytes()
            .iter()
            .enumerate()
            .map(|(i, b)| match b {
                b'A' | b'T' | b'C' | b'G' => *b as char,
                _ => ['A', 'T', 'C', 'G'][i % 4],
            })
            .collect();
        let labels = hmm.label(&seq);
        let purity = calc_purity(&seq.as_bytes(), &hmm, &motifs, &labels);
        let labels = hmm.label_motifs(&labels);
        // Remove labels corresponding to the skip state
        let labels = labels
            .into_iter()
            .filter(|rec| rec.motif_index < motifs.len())
            .collect_vec();
        let motif_counts = count_motifs(&locus.motifs, &labels);
        let labels = collapse_labels(labels);
        // TODO: Consider using empty labels instead of None
        let labels = if !labels.is_empty() {
            Some(labels)
        } else {
            None
        };

        annotations.push(Annotation {
            labels,
            motif_counts,
            purity,
        });
    }

    annotations
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
    let mut mes = Vec::new();
    mes.reserve(motifs.len() + 1);
    let mut ms = rs + 1;
    for motif in motifs {
        let num_motif_states = 3 * motif.len() + 1;
        let me = ms + num_motif_states - 1;

        hmm.set_ems(ms, vec![0.00, 0.00, 0.00, 0.00, 0.00]);
        hmm.set_trans(ms, vec![rs, me], vec![rs_to_ms, 1.0 - me_to_re]);

        define_motif_block(&mut hmm, ms, &motif);

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

fn define_motif_block(hmm: &mut Hmm, ms: usize, motif: &Vec<u8>) {
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

    // Define insersion states
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

fn get_match_emissions(char: u8) -> Vec<f64> {
    match char {
        b'A' => vec![0.00, 0.90, 0.03, 0.03, 0.03],
        b'T' => vec![0.00, 0.03, 0.90, 0.03, 0.03],
        b'C' => vec![0.00, 0.03, 0.03, 0.90, 0.03],
        b'G' => vec![0.00, 0.03, 0.03, 0.03, 0.90],
        _ => panic!("Enountered unknown base {char}"),
    }
}

#[derive(Debug, PartialEq)]
pub enum HmmEvent {
    Match,
    Mismatch,
    Ins,
    Del,
    Trans, // Silent states that don't encode alignment operation
    Skip,  // Skip state that matches bases outside of any motif run
}

pub fn get_events(hmm: &Hmm, motifs: &[Vec<u8>], states: &[usize], query: &[u8]) -> Vec<HmmEvent> {
    let mut state_to_hmm_motif = vec![-1; hmm.num_states];
    for hmm_motif in &hmm.motifs {
        // Skip two terminal states MS and ME
        for state in hmm_motif.start_state + 1..hmm_motif.end_state {
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
    for state in states {
        let motif_index = state_to_hmm_motif[*state];
        let event = if motif_index == -1 {
            HmmEvent::Trans
        } else if motif_index as usize + 1 == hmm.motifs.len() {
            HmmEvent::Skip
        } else {
            let motif = &hmm.motifs[motif_index as usize];
            let offset = state - motif.start_state - 1;
            let motif_len = motifs[motif.motif_index].len();
            match offset.div_euclid(motif_len) {
                0 => {
                    let base = query[base_index];
                    if base == get_base_match(&hmm, *state) {
                        HmmEvent::Match
                    } else {
                        HmmEvent::Mismatch
                    }
                }
                1 => HmmEvent::Ins,
                2 => HmmEvent::Del,
                _ => handle_error_and_exit(format!("Event decoding error")),
            }
        };

        if base_consumers.contains(&event) {
            base_index += 1;
        }
        events.push(event);
    }

    events
}

fn get_base_match(hmm: &Hmm, state: usize) -> u8 {
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

fn calc_purity(query: &[u8], hmm: &Hmm, motifs: &[Vec<u8>], states: &[usize]) -> f64 {
    let events = get_events(hmm, motifs, states, query);
    let num_matches = events
        .iter()
        .map(|e| match *e {
            HmmEvent::Del | HmmEvent::Ins | HmmEvent::Mismatch | HmmEvent::Skip => -1,
            HmmEvent::Match => 1,
            HmmEvent::Trans => 0,
        })
        .sum::<i32>();
    num_matches as f64 / query.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn summarize(spans: &Vec<Span>) -> Vec<(usize, usize, usize)> {
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

    #[test]
    fn calculate_purity_of_perfect_repeats() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "CAGCAGCAGCCGCCGCCGCCG";
        let states = hmm.label(&query);
        assert_eq!(calc_purity(&query.as_bytes(), &hmm, &motifs, &states), 1.0);
    }

    #[test]
    fn calculate_purity_of_imperfect_repeats() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "CAGCGCAGCCGCCGCCGGG";
        let states = hmm.label(&query);

        let purity = calc_purity(&query.as_bytes(), &hmm, &motifs, &states);
        assert_eq!(purity, 16.0 / 19.0);
    }

    #[test]
    fn calculate_purity_of_repeats_with_skip_states() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "CAGCAGCAGTTTTTTTTCCGCCGCCG";
        let states = hmm.label(&query);

        let purity = calc_purity(&query.as_bytes(), &hmm, &motifs, &states);
        assert_eq!(purity, 10.0 / 26.0);
    }
}
