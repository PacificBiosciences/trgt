use super::align::{AlignOp, AlignSeg, MotifBound};
use super::locus::Locus;
use crate::hmm::{build_hmm, get_events, HmmEvent};
use crate::trvz::align::SegType;
use itertools::Itertools;

/// Aligns a given allele to a perfect repeat as specified by the locus definition
pub fn align_consensus(locus: &Locus, consensus: &str) -> (Vec<AlignSeg>, Vec<MotifBound>) {
    let mut align = vec![AlignSeg {
        width: locus.left_flank.len(),
        op: AlignOp::Match,
        seg_type: SegType::LeftFlank,
    }];

    let motifs = locus
        .motifs
        .iter()
        .map(|m| m.as_bytes().to_vec())
        .collect_vec();
    let query = &consensus[locus.left_flank.len()..consensus.len() - locus.right_flank.len()];
    let (query_align, motif_bounds) = align_motifs(&motifs, query);
    let motif_bounds = motif_bounds
        .into_iter()
        .map(|bound| MotifBound {
            start: locus.left_flank.len() + bound.start,
            end: locus.left_flank.len() + bound.end,
            motif_index: bound.motif_index,
        })
        .collect_vec();
    align.extend(query_align);
    align.push(AlignSeg {
        width: locus.right_flank.len(),
        op: AlignOp::Match,
        seg_type: SegType::RightFlank,
    });
    (align, motif_bounds)
}

/// Aligns the given sequence to a perfect repeat composed of the given set of motifs
fn align_motifs(motifs: &[Vec<u8>], seq: &str) -> (Vec<AlignSeg>, Vec<MotifBound>) {
    if seq.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let hmm = build_hmm(motifs);
    let states = hmm.label(seq);
    let motif_spans = hmm.label_motifs(&states);
    let mut motif_by_base = vec![motifs.len(); seq.len()];
    for span in motif_spans {
        motif_by_base[span.start..span.end].fill(span.motif_index);
    }

    let events = get_events(&hmm, motifs, &states, seq.as_bytes());
    let mut align = Vec::new();
    let mut bound_events = Vec::new();
    let mut base_pos = 0;

    for (event, group) in &events.iter().chunk_by(|e| *e) {
        let seg_type = if base_pos < motif_by_base.len() {
            SegType::Tr(motif_by_base[base_pos])
        } else {
            assert_eq!(base_pos, motif_by_base.len());
            SegType::Tr(motif_by_base[base_pos.saturating_sub(1)])
        };

        let width = group.count();
        match event {
            HmmEvent::Trans => (),
            HmmEvent::MotifStart | HmmEvent::MotifEnd => {
                if let SegType::Tr(index) = seg_type {
                    bound_events.push((event, base_pos, index));
                }
            }
            HmmEvent::Del => align.push(AlignSeg {
                width: 0,
                op: AlignOp::Ins,
                seg_type,
            }),
            HmmEvent::Ins => align.push(AlignSeg {
                width,
                op: AlignOp::Del,
                seg_type,
            }),
            HmmEvent::Match => align.push(AlignSeg {
                width,
                op: AlignOp::Match,
                seg_type,
            }),
            HmmEvent::Mismatch => align.push(AlignSeg {
                width,
                op: AlignOp::Subst,
                seg_type,
            }),
            HmmEvent::Skip => {
                assert_eq!(seg_type, SegType::Tr(motifs.len()));
                align.push(AlignSeg {
                    width,
                    op: AlignOp::Match,
                    seg_type,
                })
            }
        }

        base_pos += match event {
            HmmEvent::Match | HmmEvent::Mismatch | HmmEvent::Ins | HmmEvent::Skip => width,
            HmmEvent::MotifStart | HmmEvent::MotifEnd | HmmEvent::Del | HmmEvent::Trans => 0,
        };
    }

    assert_eq!(base_pos, seq.len());

    let mut merged_align = Vec::new();
    let iter = align
        .into_iter()
        .chunk_by(|seg| (seg.op.clone(), seg.seg_type));
    for ((op, segment), group) in &iter {
        let width: usize = group.into_iter().map(|s| s.width).sum();
        merged_align.push(AlignSeg {
            width,
            op,
            seg_type: segment,
        });
    }

    let motif_bounds = bound_events
        .chunks(2)
        .map(|rec| {
            let (start, end) = rec
                .iter()
                .collect_tuple()
                .expect("Unexpected motif segmentation");
            MotifBound {
                start: start.1,
                end: end.1,
                motif_index: start.2,
            }
        })
        .collect();

    (merged_align, motif_bounds)
}
