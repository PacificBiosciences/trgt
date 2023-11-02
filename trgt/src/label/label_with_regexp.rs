use super::guess_motif_counts::guess_motif_counts;
use super::kmer_filter::KmerFilter;
use super::refine_motif_counts::refine_motif_counts;
use super::struc::{decode_regexp, Motif, Mult};
use super::{Annotation, Span, Spans};
use crate::locus::Locus;
use bio::alignment::pairwise::*;
use bio::alignment::{Alignment, AlignmentOperation as AlignOp};
use itertools::Itertools;

pub fn label_with_regexp(locus: &Locus, seqs: &Vec<String>) -> Vec<Annotation> {
    let motifs = decode_regexp(&locus.struc);
    let filter = KmerFilter::new(&motifs);
    let mut annotations = Vec::new();

    for seq in seqs {
        annotations.push(if filter.is_passing(seq) {
            let counts = guess_motif_counts(&motifs, seq);
            let (counts, purity) = refine_motif_counts(&motifs, seq, &counts);
            let spans = get_spans(seq, &motifs, &counts);
            let labels = relabel_motifs(locus, &motifs, spans);
            let motif_counts = count_motifs(&locus.motifs, &labels);
            let labels = Some(labels);
            Annotation {
                labels,
                motif_counts,
                purity,
            }
        } else {
            Annotation {
                labels: None,
                motif_counts: vec![0; locus.motifs.len()],
                purity: 0.0,
            }
        });
    }

    annotations
}

fn relabel_motifs(locus: &Locus, motifs: &[Motif], spans: Vec<Span>) -> Vec<Span> {
    spans
        .into_iter()
        .filter(|s| motifs[s.motif_index].mult == Mult::Many)
        .map(|s| Span {
            motif_index: locus
                .motifs
                .iter()
                .position(|m| m == &motifs[s.motif_index].seq)
                .unwrap(),
            start: s.start,
            end: s.end,
        })
        .collect_vec()
}

fn get_spans(seq: &str, motifs: &[Motif], motif_counts: &[usize]) -> Spans {
    let reference_spans = get_reference_spans(motifs, motif_counts);
    let reference = make_repeat(motifs, motif_counts);
    let align = align(seq, &reference);
    summarize(seq, &align, reference_spans)
}

fn align(query: &str, reference: &str) -> Alignment {
    let mut aligner = Aligner::with_capacity(query.len(), reference.len(), -5, -1, |a, b| {
        if a == b {
            1i32
        } else {
            -1i32
        }
    });

    aligner.global(query.as_bytes(), reference.as_bytes())
}

fn make_repeat(motifs: &[Motif], motif_counts: &[usize]) -> String {
    let mut repeat = String::new();
    for (index, motif) in motifs.iter().enumerate() {
        if motif.mult == Mult::Many {
            repeat += &motif.seq.repeat(motif_counts[index]);
        } else {
            repeat += &motif.seq;
        }
    }

    repeat
}

fn get_reference_spans(motifs: &[Motif], motif_counts: &[usize]) -> Spans {
    let mut spans = Vec::new();
    let mut last_segment_end = 0;
    for (motif_index, motif) in motifs.iter().enumerate() {
        let mult = if motif.mult == Mult::Once {
            1
        } else {
            motif_counts[motif_index]
        };

        let start = last_segment_end;
        let end = start + motif.seq.len() * mult;
        spans.push(Span {
            motif_index,
            start,
            end,
        });
        last_segment_end = end;
    }

    spans
}

fn summarize(query: &str, align: &Alignment, ref_spans: Spans) -> Spans {
    assert!(align.xstart == 0 && align.ystart == 0);
    let mut x_pos = 0; // Query pos
    let mut y_pos = 0; // Reference pos
    let mut query_motif_indexes = vec![-1; query.len()];

    for op in &align.operations {
        query_motif_indexes[x_pos] = get_motif_index(y_pos, &ref_spans);

        x_pos += match *op {
            AlignOp::Match => 1,
            AlignOp::Subst => 1,
            AlignOp::Del => 0,
            AlignOp::Ins => 1,
            _ => panic!("Unhandled operation {:?}", *op),
        };

        y_pos += match *op {
            AlignOp::Match => 1,
            AlignOp::Subst => 1,
            AlignOp::Del => 1,
            AlignOp::Ins => 0,
            _ => panic!("Unhandled operation {:?}", *op),
        };

        if x_pos == query.len() {
            break;
        }
    }

    let mut query_spans = Vec::new();
    let mut last_segment_end = 0;
    for (index, group) in &query_motif_indexes.iter().group_by(|i| *i) {
        if *index == -1 {
            last_segment_end += group.count();
            continue;
        }

        let start = last_segment_end;
        let end = start + group.count();
        query_spans.push(Span {
            motif_index: *index as usize,
            start,
            end,
        });
        last_segment_end = end;
    }

    // Add placeholder spans for missing motifs
    let mut complete_query_spans = Vec::new();
    let last_motif_index = ref_spans.last().unwrap().motif_index;
    let mut query_span_index = 0;
    for motif_index in 0..last_motif_index + 1 {
        if query_span_index < query_spans.len()
            && motif_index == query_spans[query_span_index].motif_index
        {
            complete_query_spans.push(query_spans[query_span_index].clone());
            query_span_index += 1;
        } else {
            let start = if complete_query_spans.is_empty() {
                0
            } else {
                complete_query_spans.last().unwrap().end
            };

            complete_query_spans.push(Span {
                motif_index,
                start,
                end: start,
            });
        }
    }
    complete_query_spans
}

fn get_motif_index(pos: usize, spans: &Vec<Span>) -> i32 {
    for span in spans {
        if span.start <= pos && pos < span.end {
            return span.motif_index as i32;
        }
    }

    -1
}

fn count_motifs(motifs: &Vec<String>, labels: &Spans) -> Vec<usize> {
    let mut motif_counts = vec![0; motifs.len()];

    for span in labels {
        motif_counts[span.motif_index] += span.len();
    }

    motif_counts = motif_counts
        .into_iter()
        .enumerate()
        .map(|(i, c)| (c as f64 / motifs[i].len() as f64).round() as usize)
        .collect_vec();

    motif_counts
}
