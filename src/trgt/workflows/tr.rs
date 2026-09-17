use super::{Allele, Genotype, LocusResult};
use crate::hmm::{
    build_hmm, calc_purity, collapse_labels, count_motifs, remove_imperfect_motifs,
    replace_invalid_bases, Annotation,
};
use crate::trgt::{
    genotype::{find_tr_spans, genotype_cluster, genotype_flank, genotype_size, Gt},
    locus::Locus,
    reads::{AlleleAssign, LocusRead, SpanningRead},
};
use crate::utils::{fnv1a64, Genotyper, Ploidy, Result};
use itertools::Itertools;
use rand::{seq::SliceRandom, SeedableRng};
use std::collections::BTreeMap;
use std::vec;

pub struct Params {
    pub min_flank_id_frac: f64,
    pub min_read_qual: f32,
    pub search_flank_len: usize,
    pub max_depth: usize,
    pub skip_phase_annotation: bool,
}

pub fn analyze(mut locus: Locus, params: &Params) -> (Locus, Result<LocusResult>) {
    if locus.ploidy == Ploidy::Zero {
        return (locus, Ok(LocusResult::empty()));
    }

    let spanning_reads = get_spanning_reads(&mut locus, params);

    const MIN_RQ_FOR_PURITY: f32 = 0.9;
    let mut spanning_reads = if params.min_read_qual < MIN_RQ_FOR_PURITY {
        let read_len = spanning_reads.len();
        let ret = filter_impure_trs(&locus, spanning_reads, MIN_RQ_FOR_PURITY);
        if ret.len() < read_len {
            log::warn!(
                "{}: Filtered out {} impure reads",
                locus.id,
                read_len - ret.len()
            );
        }
        ret
    } else {
        spanning_reads
    };

    log::debug!(
        "{}: Found {} spanning reads",
        locus.id,
        spanning_reads.len()
    );

    if spanning_reads.is_empty() {
        return (locus, Ok(LocusResult::empty()));
    }

    // Sort and (possibly) downsample
    spanning_reads.sort_by(|a, b| {
        a.tr_len()
            .cmp(&b.tr_len())
            .then_with(|| a.read.id.cmp(&b.read.id))
    });

    if spanning_reads.len() > params.max_depth {
        uniform_downsample(&mut spanning_reads, params.max_depth);
        log::debug!(
            "{}: downsampled to {} reads",
            locus.id,
            spanning_reads.len()
        );
    }

    deterministic_shuffle(&mut spanning_reads, &locus.id);
    normalize_by_phase_set(&mut spanning_reads);

    let trs: Vec<&[u8]> = spanning_reads.iter().map(|sr| sr.tr_slice()).collect();
    let (mut gt, mut allele_seqs, mut classification) = match locus.genotyper {
        Genotyper::Size => genotype_size::genotype(locus.ploidy, &trs),
        Genotyper::Cluster => genotype_cluster::genotype(locus.ploidy, &trs),
    };

    // Attempt flank re-genotyping only if alleles have similar length
    if gt.len() == 2 && gt[0].size.abs_diff(gt[1].size) <= 10 {
        let read_refs: Vec<&LocusRead> = spanning_reads.iter().map(|sr| &sr.read).collect();
        let snp_result = genotype_flank::genotype(&read_refs, &trs, &locus.region);
        if let Some((snp_gt, snp_alleles, snp_assignment)) = snp_result {
            (gt, allele_seqs, classification) = (snp_gt, snp_alleles, snp_assignment);
        }
    }

    for (sr, cls) in spanning_reads.iter_mut().zip(classification.iter()) {
        sr.assign = match *cls {
            0 => AlleleAssign::A0,
            1 => AlleleAssign::A1,
            _ => AlleleAssign::None,
        };
    }

    let annotations = label_with_hmm(&locus, &allele_seqs);

    let spanning_by_hap = [
        spanning_reads
            .iter()
            .filter(|sr| matches!(sr.assign, AlleleAssign::A0 | AlleleAssign::Both))
            .count(),
        spanning_reads
            .iter()
            .filter(|sr| matches!(sr.assign, AlleleAssign::A1 | AlleleAssign::Both))
            .count(),
    ];

    let meth_by_hap = get_meth_by_hap(&gt, &spanning_reads);

    let mut genotype = Genotype::new();
    for (i, (g, (seq, annotation))) in gt
        .into_iter()
        .zip(allele_seqs.into_iter().zip(annotations.into_iter()))
        .enumerate()
    {
        genotype.push(Allele {
            seq,
            annotation,
            ci: g.ci,
            num_spanning: spanning_by_hap[i],
            meth: meth_by_hap[i],
        });
    }

    // Put reference allele first
    if genotype.len() != 1 && genotype[0].seq != locus.tr && genotype[1].seq == locus.tr {
        genotype.swap(0, 1);
        for sr in &mut spanning_reads {
            sr.assign = match sr.assign {
                AlleleAssign::A0 => AlleleAssign::A1,
                AlleleAssign::A1 => AlleleAssign::A0,
                x => x,
            };
        }
    }

    let phase_result = if params.skip_phase_annotation {
        None
    } else {
        phase_by_hp_ps(&spanning_reads)
    };

    (
        locus,
        Ok(LocusResult {
            genotype,
            spanning_reads,
            phase_result,
        }),
    )
}

fn get_spanning_reads(locus: &mut Locus, params: &Params) -> Vec<SpanningRead> {
    let reads = std::mem::take(&mut locus.reads);
    let tr_spans = find_tr_spans(&locus.left_flank, &locus.right_flank, &reads, params);
    reads
        .into_iter()
        .zip(tr_spans)
        .filter_map(|(read, span)| match span {
            Some((s, e))
                if s >= params.search_flank_len
                    && read.bases.len() - e >= params.search_flank_len =>
            {
                Some(SpanningRead {
                    read,
                    span: (s, e),
                    assign: AlleleAssign::None,
                })
            }
            _ => None,
        })
        .collect()
}

fn uniform_downsample(spanning_reads: &mut Vec<SpanningRead>, output_len: usize) {
    let n = spanning_reads.len() as f64;
    let mut fast: f64 = 0.0;
    let step = n / (output_len as f64);
    for i in 0..output_len {
        let ind = fast.floor() as usize;
        if ind != i {
            spanning_reads.swap(i, ind);
        }
        fast += step;
    }
    spanning_reads.truncate(output_len);
}

#[inline]
fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
}

fn get_meth_by_hap(gt: &Gt, spanning_reads: &[SpanningRead]) -> Vec<Option<f64>> {
    let mut meths_1 = Vec::new();
    let mut meths_2 = Vec::new();

    for sr in spanning_reads {
        if let Some(level) = get_tr_meth(sr) {
            match assign_read(gt, sr.tr_len()) {
                Assignment::First => meths_1.push(level),
                Assignment::Second => meths_2.push(level),
                Assignment::Both => {
                    meths_1.push(level);
                    meths_2.push(level);
                }
                Assignment::None => {}
            };
        }
    }

    if gt.len() == 2 {
        vec![mean(&meths_1), mean(&meths_2)]
    } else {
        vec![mean(&meths_1)]
    }
}

enum Assignment {
    First,
    Second,
    Both,
    None,
}

fn assign_read(gt: &Gt, tr_len: usize) -> Assignment {
    if gt.len() == 1 {
        return Assignment::First;
    }

    let hap1_len = gt[0].size;
    let hap2_len = gt[1].size;

    let spans_1 = gt[0].ci.0 <= tr_len && tr_len <= gt[0].ci.1;
    let spans_2 = gt[1].ci.0 <= tr_len && tr_len <= gt[1].ci.1;

    let dist_1 = tr_len.abs_diff(gt[0].size);
    let dist_2 = tr_len.abs_diff(gt[1].size);

    if dist_1 < dist_2 && spans_1 {
        return Assignment::First;
    }

    if dist_2 < dist_1 && spans_2 {
        return Assignment::Second;
    }

    if hap1_len == hap2_len && spans_1 {
        return Assignment::Both;
    }

    Assignment::None
}

fn get_tr_meth(sr: &SpanningRead) -> Option<f64> {
    let meth = match &sr.read.meth {
        Some(m) if !m.is_empty() => m,
        _ => return None,
    };

    let mut total_meth = 0.0;
    let mut cpg_count = 0;
    let mut cpg_index = 0;
    for (pos, dinuc) in sr.read.bases.windows(2).enumerate() {
        if dinuc == b"CG" {
            if sr.span.0 <= pos && pos < sr.span.1 {
                cpg_count += 1;
                total_meth += match meth.get(cpg_index) {
                    Some(value) => *value as f64 / 255.0,
                    None => {
                        log::error!("Read {} has malformed methylation profile", sr.read.id);
                        std::process::exit(1);
                    }
                };
            }

            cpg_index += 1;
        }
    }

    if cpg_count != 0 {
        let mean_meth = total_meth / cpg_count as f64;
        Some(mean_meth)
    } else {
        None
    }
}

fn filter_impure_trs(
    locus: &Locus,
    spanning_reads: Vec<SpanningRead>,
    rq_cutoff: f32,
) -> Vec<SpanningRead> {
    const PURITY_CUTOFF: f64 = 0.9;
    let n = spanning_reads.len();
    if n == 0 {
        return spanning_reads;
    }

    let low_rq_indices: Vec<usize> = spanning_reads
        .iter()
        .enumerate()
        .filter_map(|(i, sr)| {
            if sr.read.read_qual.is_none_or(|rq| rq < rq_cutoff) {
                Some(i)
            } else {
                None
            }
        })
        .collect();

    if low_rq_indices.is_empty() {
        return spanning_reads;
    }

    let motifs: Vec<_> = locus
        .motifs
        .iter()
        .map(|m| replace_invalid_bases(m, b"ATCGN"))
        .collect();
    let hmm = build_hmm(&motifs);

    let mut candidates: Vec<(usize, f64)> = low_rq_indices
        .into_iter()
        .filter_map(|i| {
            let seq = replace_invalid_bases(spanning_reads[i].tr_slice(), b"ATCG");
            let labels = hmm.label(&seq);
            let p = calc_purity(&seq, &hmm, &motifs, &labels);
            if p < PURITY_CUTOFF {
                Some((i, p))
            } else {
                None
            }
        })
        .collect();

    let max_filter = std::cmp::max(1_usize, (0.1 * n as f64).round() as usize);
    let num_to_drop = std::cmp::min(max_filter, candidates.len());
    let mut drop_mask = vec![false; n];
    if num_to_drop > 0 {
        if candidates.len() <= num_to_drop {
            for (i, _) in candidates {
                drop_mask[i] = true;
            }
        } else {
            candidates.select_nth_unstable_by(num_to_drop - 1, |a, b| {
                a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0))
            });
            for (i, _) in candidates.iter().take(num_to_drop) {
                drop_mask[*i] = true;
            }
        }
    }

    let mut new_spanning_reads = Vec::with_capacity(n - num_to_drop);
    for (i, sr) in spanning_reads.into_iter().enumerate() {
        if !drop_mask[i] {
            new_spanning_reads.push(sr);
        }
    }

    new_spanning_reads
}

fn label_with_hmm(locus: &Locus, seqs: &[Vec<u8>]) -> Vec<Annotation> {
    let motifs = locus
        .motifs
        .iter()
        .map(|m| replace_invalid_bases(m, b"ATCGN"))
        .collect_vec();
    let hmm = build_hmm(&motifs);

    let mut annotations = Vec::new();
    for seq in seqs {
        let seq = replace_invalid_bases(seq, b"ATCG");
        let labels = hmm.label(&seq);
        let purity = calc_purity(&seq, &hmm, &motifs, &labels);
        let labels = remove_imperfect_motifs(&hmm, &motifs, &labels, &seq, 6);
        let labels = hmm.label_motifs(&labels);
        // Remove labels corresponding to the skip state
        let labels = labels
            .into_iter()
            .filter(|rec| rec.motif_index < motifs.len())
            .collect_vec();
        let motif_counts = count_motifs(&locus.motifs, &labels);
        let labels = collapse_labels(labels);
        // TODO: Consider using empty labels instead of Nones
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

#[inline]
pub fn deterministic_shuffle(reads: &mut [SpanningRead], locus_id: &str) {
    if reads.len() <= 1 {
        return;
    }
    let mut rng = rand_chacha::ChaCha12Rng::seed_from_u64(fnv1a64(locus_id.as_bytes()));
    reads.shuffle(&mut rng)
}

fn normalize_by_phase_set(reads: &mut [SpanningRead]) {
    let mut ps_counts: BTreeMap<u32, usize> = BTreeMap::new();
    for sr in reads.iter().filter(|sr| sr.read.ps_tag.is_some()) {
        *ps_counts.entry(sr.read.ps_tag.unwrap()).or_insert(0) += 1;
    }

    let Some((&dominant_ps, _)) = ps_counts.iter().max_by(|(ps_a, count_a), (ps_b, count_b)| {
        count_a.cmp(count_b).then_with(|| ps_b.cmp(ps_a))
    }) else {
        return;
    };

    for sr in reads {
        if sr.read.ps_tag != Some(dominant_ps) {
            sr.read.ps_tag = None;
            sr.read.hp_tag = None;
        }
    }
}

#[derive(Debug)]
pub struct PhaseResult {
    /// Phased genotype order to emit in VCF GT, e.g. (0,1) means "0|1".
    pub phased_gt_digits: (u8, u8),
    pub ps: Option<u32>,
}

pub fn phase_by_hp_ps(reads: &[SpanningRead]) -> Option<PhaseResult> {
    const HP_THRESHOLD: f64 = 0.7;
    const CONFLICT_THRESHOLD: f64 = 0.1;

    if reads.is_empty() {
        return None;
    }

    let mut chosen_ps = None;
    for r in reads.iter() {
        if let Some(ps) = r.read.ps_tag {
            match chosen_ps {
                Some(existing) if existing != ps => {
                    return None;
                }
                None => chosen_ps = Some(ps),
                _ => {}
            }
        }
    }

    let mut c = [[0usize; 2]; 2];
    let mut n_reads = 0usize;
    let mut n_reads_with_hp = 0usize;
    let mut seen_hp = [false; 2];
    for r in reads.iter() {
        n_reads += 1;
        let Some(hp) = r.read.hp_tag else { continue };
        debug_assert!(hp == 1 || hp == 2);
        let hp_idx = (hp - 1) as usize;
        n_reads_with_hp += 1;
        let al = match r.assign {
            AlleleAssign::A0 => 0,
            AlleleAssign::A1 => 1,
            _ => continue,
        };
        seen_hp[hp_idx] = true;
        c[hp_idx][al] += 1;
    }

    if !seen_hp.iter().all(|seen| *seen) {
        return None;
    }

    let hp_fraction = n_reads_with_hp as f64 / n_reads as f64;
    if hp_fraction < HP_THRESHOLD {
        return None;
    }

    // A: HP1->allele0, HP2->allele1  (gives GT "0|1")
    // B: HP1->allele1, HP2->allele0  (gives GT "1|0")
    let agree_a = c[0][0] + c[1][1]; // 0|1
    let agree_b = c[0][1] + c[1][0]; // 1|0

    let conflict_fraction = 1.0 - (agree_a.max(agree_b) as f64 / n_reads_with_hp as f64);
    if conflict_fraction > CONFLICT_THRESHOLD {
        return None;
    }

    // Orientation that maximises agreement.
    let phased_gt_digits = if agree_a >= agree_b {
        (0u8, 1u8)
    } else {
        (1u8, 0u8)
    };

    Some(PhaseResult {
        phased_gt_digits,
        ps: chosen_ps,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::trgt::reads::HiFiRead;

    fn build_read(id: usize, hp: Option<u8>, ps: Option<u32>) -> SpanningRead {
        SpanningRead {
            read: HiFiRead {
                id: format!("read-{id}"),
                is_reverse: false,
                bases: vec![b'A'; 10],
                quals: vec![30; 10],
                meth: None,
                read_qual: None,
                mismatch_positions: None,
                cigar: None,
                hp_tag: hp,
                ps_tag: ps,
                mapq: 60,
                ref_start: 0,
                ref_end: 10,
            },
            span: (0, 2),
            assign: AlleleAssign::None,
        }
    }

    #[test]
    fn normalize_reads_keeps_dominant_ps() {
        let mut reads = vec![
            build_read(1, Some(1), Some(10)),
            build_read(2, Some(1), Some(10)),
            build_read(3, Some(2), Some(20)),
            build_read(4, None, None),
        ];

        normalize_by_phase_set(&mut reads);

        assert_eq!(reads[0].read.ps_tag, Some(10));
        assert_eq!(reads[0].read.hp_tag, Some(1));
        assert_eq!(reads[1].read.ps_tag, Some(10));
        assert_eq!(reads[1].read.hp_tag, Some(1));
        assert_eq!(reads[2].read.ps_tag, None);
        assert_eq!(reads[2].read.hp_tag, None);
        assert_eq!(reads[3].read.ps_tag, None);
        assert_eq!(reads[3].read.hp_tag, None);
    }

    #[test]
    fn normalize_reads_prefers_smaller_ps_on_ties() {
        let mut reads = vec![
            build_read(1, Some(1), Some(10)),
            build_read(2, Some(1), Some(20)),
            build_read(3, Some(2), Some(20)),
            build_read(4, Some(2), Some(10)),
        ];

        normalize_by_phase_set(&mut reads);

        assert_eq!(reads[0].read.ps_tag, Some(10));
        assert_eq!(reads[1].read.ps_tag, None);
        assert_eq!(reads[2].read.ps_tag, None);
        assert_eq!(reads[3].read.ps_tag, Some(10));
    }

    #[test]
    fn normalize_reads_no_ps_observed() {
        let mut reads = vec![
            build_read(1, Some(1), None),
            build_read(2, None, None),
            build_read(3, None, None),
        ];

        normalize_by_phase_set(&mut reads);

        assert_eq!(reads[0].read.ps_tag, None);
        assert_eq!(reads[0].read.hp_tag, Some(1));
        assert_eq!(reads[1].read.ps_tag, None);
        assert_eq!(reads[1].read.hp_tag, None);
        assert_eq!(reads[2].read.ps_tag, None);
        assert_eq!(reads[2].read.hp_tag, None);
    }

    #[test]
    fn phase_by_hp_ps_rejects_multiple_ps() {
        let mut reads = vec![
            build_read(1, Some(1), Some(10)),
            build_read(2, Some(2), Some(20)),
        ];
        for sr in &mut reads {
            sr.assign = AlleleAssign::A0;
        }

        assert!(phase_by_hp_ps(&reads).is_none());
    }

    #[test]
    fn phase_by_hp_ps_requires_two_hps() {
        let mut reads = vec![
            build_read(1, Some(1), Some(10)),
            build_read(2, Some(1), Some(10)),
        ];
        reads[0].assign = AlleleAssign::A0;
        reads[1].assign = AlleleAssign::A1;

        assert!(phase_by_hp_ps(&reads).is_none());
    }

    #[test]
    fn phase_by_hp_ps_counts_hp_reads_when_assigned() {
        let mut reads = vec![
            build_read(1, Some(1), Some(10)),
            build_read(2, Some(1), Some(10)),
            build_read(3, Some(2), Some(10)),
            build_read(4, Some(2), Some(10)),
        ];
        let assigns = [
            AlleleAssign::A0,
            AlleleAssign::A0,
            AlleleAssign::A1,
            AlleleAssign::A1,
        ];
        for (sr, assign) in reads.iter_mut().zip(assigns) {
            sr.assign = assign;
        }

        let pr = phase_by_hp_ps(&reads).expect("phase result should be produced");
        assert_eq!(pr.phased_gt_digits, (0, 1));
        assert_eq!(pr.ps, Some(10));
    }
}
