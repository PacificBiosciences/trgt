use crate::input::Spans;
use crate::locus::{BaseLabel, Locus};
use crate::struc::{self, RegionLabel};
use crate::struc::{Mult, Seq};
use bio::alignment::{pairwise::*, AlignmentOperation};
use itertools::Itertools;

pub fn label_with_motifs(
    locus: &Locus,
    spans_by_allele: &[Option<Spans>],
    alleles: &[String],
) -> Vec<Vec<BaseLabel>> {
    let tr_ref_spans = get_tr_ref_spans(locus, spans_by_allele);
    let ref_alleles = get_ref_alleles(locus, alleles, spans_by_allele);

    let mut base_labels_by_allele = Vec::new();
    for (index, ref_seq) in ref_alleles.iter().enumerate() {
        let allele = &alleles[index];
        let tr_spans = &tr_ref_spans[index];
        base_labels_by_allele.push(label_bases(locus, allele, ref_seq, tr_spans));
    }
    base_labels_by_allele
}

fn label_bases(
    locus: &Locus,
    allele: &str,
    ref_seq: &str,
    ref_tr_spans: &Option<Vec<RegionLabel>>,
) -> Vec<BaseLabel> {
    let mut aligner = Aligner::with_capacity(allele.len(), ref_seq.len(), -5, -1, |a, b| {
        if a == b {
            1i32
        } else {
            -1i32
        }
    });

    let allele = &allele[locus.left_flank.len()..allele.len() - locus.right_flank.len()];
    let ref_seq = &ref_seq[locus.left_flank.len()..ref_seq.len() - locus.right_flank.len()];

    let align = aligner.global(allele.as_bytes(), ref_seq.as_bytes());

    let mut base_labels = vec![BaseLabel::Match; locus.left_flank.len()];
    let mut ref_pos = locus.left_flank.len();
    for (op_index, op) in align.operations.iter().enumerate() {
        if let Some(index) = get_span_index(ref_tr_spans, ref_pos) {
            // TODO: Consolidate
            if let RegionLabel::Tr(start, _end, motif) = &ref_tr_spans.as_ref().unwrap()[index] {
                if (ref_pos - start) % motif.len() == 0
                    && (op_index == 0 || align.operations[op_index - 1] != AlignmentOperation::Ins)
                {
                    base_labels.push(BaseLabel::MotifBound);
                }
            } else {
                panic!();
            }
        }

        base_labels.push(match op {
            AlignmentOperation::Del => BaseLabel::Skip,
            AlignmentOperation::Ins => BaseLabel::NoMatch,
            AlignmentOperation::Match => BaseLabel::Match,
            AlignmentOperation::Subst => BaseLabel::Mismatch,
            _ => panic!("Unable to handle {:?}", op),
        });

        // allele_pos += match op {
        //     AlignmentOperation::Del => 0,
        //    AlignmentOperation::Ins => 1,
        //    AlignmentOperation::Match => 1,
        //    AlignmentOperation::Subst => 1,
        //    _ => panic!("Unable to handle {:?}", op),
        //};

        ref_pos += match op {
            AlignmentOperation::Del => 1,
            AlignmentOperation::Ins => 0,
            AlignmentOperation::Match => 1,
            AlignmentOperation::Subst => 1,
            _ => panic!("Unable to handle {:?}", op),
        };
    }

    base_labels.extend(vec![BaseLabel::Match; locus.right_flank.len()]);
    base_labels
}

fn get_span_index(ref_tr_spans: &Option<Vec<RegionLabel>>, pos: usize) -> Option<usize> {
    if let Some(ref_tr_spans) = ref_tr_spans {
        for (index, span) in ref_tr_spans.iter().enumerate() {
            let (start, end) = match span {
                RegionLabel::Tr(start, end, _) => (start, end),
                _ => panic!(),
            };

            // End is inclusive to account for the end of the last motif
            if *start <= pos && pos <= *end {
                return Some(index);
            }
        }
    }

    None
}

fn get_ref_alleles(
    locus: &Locus,
    allele_seqs: &[String],
    spans_by_allele: &[Option<Spans>],
) -> Vec<String> {
    let mut ref_seqs = Vec::new();
    for (index, spans) in spans_by_allele.iter().enumerate() {
        ref_seqs.push(match spans {
            Some(spans) => get_ref_seq(locus, spans),
            None => allele_seqs[index].to_string(),
        });
    }
    ref_seqs
}

fn get_ref_seq(locus: &Locus, spans: &Spans) -> String {
    let mut allele_seq = locus.left_flank.clone();

    let mut tr_index = 0;
    let struc = struc::struc(&locus.struc);

    for seq in struc {
        match seq.mult {
            Mult::Once => allele_seq += &seq.motif,
            Mult::Many => {
                let span = &spans[tr_index];
                let num_motifs = (span.end - span.start) / seq.motif.len();
                let leftover_len = (span.end - span.start) % seq.motif.len();
                allele_seq += &seq.motif.clone().repeat(num_motifs);
                allele_seq += &seq.motif[..leftover_len];
                tr_index += 1;
            }
        }
    }

    allele_seq += &locus.right_flank;
    allele_seq
}

fn get_tr_ref_spans(
    locus: &Locus,
    spans_by_allele: &[Option<Spans>],
) -> Vec<Option<Vec<RegionLabel>>> {
    let struc = struc::struc(&locus.struc);

    spans_by_allele
        .iter()
        .map(|spans| match spans {
            None => None,
            Some(spans) => Some(get_tr_spans(locus, &struc, spans)),
        })
        .collect_vec()
}

fn get_tr_spans(locus: &Locus, struc: &Vec<Seq>, spans: &Spans) -> Vec<RegionLabel> {
    let mut tr_spans = Vec::new();
    let mut ref_pos = locus.left_flank.len();
    let mut tr_index = 0;
    for seq in struc {
        if seq.mult == Mult::Many {
            let tr_len = spans[tr_index].end - spans[tr_index].start;
            tr_spans.push(RegionLabel::Tr(
                ref_pos,
                ref_pos + tr_len,
                seq.motif.clone(),
            ));
            tr_index += 1;
            ref_pos += tr_len;
        } else {
            ref_pos += seq.motif.len();
        }
    }

    tr_spans
}
