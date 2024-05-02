use super::{consensus, Gt, TrSize};
use crate::trgt::reads::HiFiRead;
use crate::utils::{align, median};
use itertools::Itertools;
use std::{cmp::Ordering, collections::BTreeMap};

type Profile = Vec<Option<bool>>;

pub fn genotype(reads: &[HiFiRead], tr_seqs: &[&str]) -> Option<(Gt, Vec<String>, Vec<i32>)> {
    let (trs_by_allele, mut allele_assignment) =
        get_trs_with_hp(reads, tr_seqs).or_else(|| get_trs_with_clustering(reads, tr_seqs))?;
    let mut gt = Gt::new();
    let mut alleles = Vec::new();

    for trs in trs_by_allele {
        let (backbone, frequency) = simple_consensus(&trs)?;
        const MIN_FREQ_TO_ALIGN: f64 = 0.5;
        let allele = if frequency < MIN_FREQ_TO_ALIGN {
            let aligns = align(&backbone, &trs);
            consensus::repair_consensus(&backbone, &trs, &aligns)
        } else {
            backbone.to_string()
        };

        let min_tr_len = trs.iter().map(|tr| tr.len()).min().unwrap();
        let max_tr_len = trs.iter().map(|tr| tr.len()).max().unwrap();

        let size = TrSize::new(allele.len(), (min_tr_len, max_tr_len));
        gt.push(size);
        alleles.push(allele);
    }

    // Smaller allele should always appear first
    if alleles[0].len() > alleles[1].len() {
        gt.swap(0, 1);
        alleles.swap(0, 1);
        allele_assignment = allele_assignment.into_iter().map(|a| (a + 1) % 2).collect();
    }

    Some((gt, alleles, allele_assignment))
}

fn get_trs_with_hp<'a>(
    reads: &[HiFiRead],
    tr_seqs: &[&'a str],
) -> Option<([Vec<&'a str>; 2], Vec<i32>)> {
    let mut allele_assignment = Vec::new();
    let mut trs_by_allele = [Vec::new(), Vec::new()];
    let mut assignment_tie_breaker: usize = 1;
    let mut num_unassigned = 0;
    for (read, tr_seq) in reads.iter().zip(tr_seqs.iter()) {
        match read.hp_tag {
            Some(1) => {
                allele_assignment.push(0_i32);
                trs_by_allele[0].push(*tr_seq);
            }
            Some(2) => {
                allele_assignment.push(1_i32);
                trs_by_allele[1].push(*tr_seq);
            }
            _ => {
                assignment_tie_breaker = (assignment_tie_breaker + 1) % 2;
                allele_assignment.push(assignment_tie_breaker as i32);
                trs_by_allele[assignment_tie_breaker].push(*tr_seq);
                num_unassigned += 1;
            }
        }
    }

    let prop_assigned = (reads.len() - num_unassigned) as f64 / reads.len() as f64;
    if !trs_by_allele[0].is_empty() && !trs_by_allele[1].is_empty() && prop_assigned >= 0.7 {
        Some((trs_by_allele, allele_assignment))
    } else {
        None
    }
}

fn get_trs_with_clustering<'a>(
    reads: &[HiFiRead],
    tr_seqs: &[&'a str],
) -> Option<([Vec<&'a str>; 2], Vec<i32>)> {
    if tr_seqs.is_empty() {
        return None;
    }

    let analysis_region = get_analysis_region(reads);

    let snvs = call_snvs(analysis_region, reads, 0.20);
    let profiles = get_profiles(reads, &snvs);
    let candidate_gts = get_candidate_gts(&profiles);

    // if there is just one genotype, it must be homozygous
    if candidate_gts.len() <= 1 {
        return None;
    }

    let gt_logliks = candidate_gts
        .iter()
        .map(|gt| get_loglik(gt, &profiles))
        .collect_vec();

    let top_gt = &candidate_gts
        .iter()
        .zip(gt_logliks.iter())
        .max_by(|(_gt1, ll1), (_gt2, ll2)| ll1.partial_cmp(ll2).unwrap())
        .unwrap()
        .0;

    let mut allele_assignment = Vec::new();
    let mut assignment_tie_breaker = 1;
    let mut trs_by_allele = [Vec::new(), Vec::new()];
    for (index, profile) in profiles.iter().enumerate() {
        let dist1 = get_dist(profile, &top_gt.0);
        let dist2 = get_dist(profile, &top_gt.1);

        match dist1.cmp(&dist2) {
            Ordering::Less => {
                allele_assignment.push(0);
                trs_by_allele[0].push(tr_seqs[index]);
            }
            Ordering::Greater => {
                allele_assignment.push(1);
                trs_by_allele[1].push(tr_seqs[index]);
            }
            Ordering::Equal => {
                assignment_tie_breaker = (assignment_tie_breaker + 1) % 2;
                allele_assignment.push(assignment_tie_breaker);
                trs_by_allele[0].push(tr_seqs[index]);
                trs_by_allele[1].push(tr_seqs[index]);
            }
        }
    }
    Some((trs_by_allele, allele_assignment))
}

fn get_dist(read: &[Option<bool>], allele: &[bool]) -> usize {
    read.iter()
        .zip(allele.iter())
        .map(|(p, h)| (p.as_ref() == Some(h)) as usize)
        .sum()
}

/// Determine consensus for the input sequences
///
/// Return the most frequent sequence and its relative frequency.
/// If multiple sequences have the same frequency,
/// return the one whose length is closest to the median.
///
fn simple_consensus(seqs: &[&str]) -> Option<(String, f64)> {
    let median_len = median(&seqs.iter().map(|s| s.len() as i32).collect_vec())? as usize;
    let mut seq_to_count = BTreeMap::new();
    for seq in seqs {
        *seq_to_count.entry(*seq).or_insert(0) += 1;
    }
    let top_group_size = *seq_to_count.values().max()?;
    let consensus = seq_to_count
        .into_iter()
        .filter(|(_, c)| *c == top_group_size)
        .map(|(s, _)| (s, s.len().abs_diff(median_len)))
        .min_by_key(|(_, delta)| *delta)
        .map(|(s, _)| s.to_string())?;
    let top_group_frequency = (top_group_size as f64) / (seqs.len() as f64);
    Some((consensus, top_group_frequency))
}

fn get_loglik(gt: &(Vec<bool>, Vec<bool>), profiles: &Vec<Profile>) -> f64 {
    let mut total_ll = 0.0;
    for profile in profiles {
        let term1 = eval_profile_given_hap(profile, &gt.0);
        let term2 = eval_profile_given_hap(profile, &gt.1);
        let term = ln_sum_exp(term1, term2) - 2.0_f64.ln();
        total_ll += term;
    }

    total_ll
}

fn eval_profile_given_hap(profile: &Profile, hap: &[bool]) -> f64 {
    const MATCH_PROB: f64 = 0.9;
    const MISMATCH_PROB: f64 = 1.0 - MATCH_PROB;
    profile
        .iter()
        .zip(hap.iter())
        .filter_map(|(p, h)| {
            p.as_ref().map(|p_val| {
                if p_val == h {
                    MATCH_PROB.ln()
                } else {
                    MISMATCH_PROB.ln()
                }
            })
        })
        .sum()
}

fn ln_sum_exp(term1: f64, term2: f64) -> f64 {
    let max_term = term1.max(term2);
    max_term + ((term1 - max_term).exp() + (term2 - max_term).exp()).ln()
}

fn get_analysis_region(reads: &[HiFiRead]) -> (i32, i32) {
    // min fraction of reads that must fully cover the left (right) flank
    const COV_READ_FRAC: f64 = 0.85;
    let skip_count = (reads.len() as f64 * (1.0 - COV_READ_FRAC)).round() as usize;

    let start = reads
        .iter()
        .map(|r| r.start_offset)
        .sorted()
        .nth_back(skip_count)
        .unwrap();

    let end = reads
        .iter()
        .map(|r| r.end_offset)
        .sorted()
        .nth(skip_count)
        .unwrap();

    (start, end)
}

fn get_candidate_gts(profiles: &[Profile]) -> Vec<(Vec<bool>, Vec<bool>)> {
    let haps = profiles
        .iter()
        .filter(|p| p.iter().all(Option::is_some))
        .sorted()
        .collect_vec();

    // min fraction of reads that qualify as putative haplotypes
    const PUTATIVE_HAP_FRAC: f64 = 0.40;
    if (haps.len() as f64 / profiles.len() as f64) < PUTATIVE_HAP_FRAC {
        return Vec::new();
    }

    let haps = haps.into_iter().dedup().collect_vec();
    haps.iter()
        .enumerate()
        .flat_map(|(i, hap1)| {
            let hap1_unwrapped = hap1.iter().filter_map(|v| *v).collect_vec();
            haps[i..].iter().map(move |hap2| {
                let hap2_unwrapped = hap2.iter().filter_map(|v| *v).collect_vec();
                (hap1_unwrapped.clone(), hap2_unwrapped)
            })
        })
        .collect()
}

fn get_profiles(reads: &[HiFiRead], snvs: &[i32]) -> Vec<Profile> {
    let profiles = reads
        .iter()
        .map(|read| match &read.mismatch_offsets {
            Some(mm_offsets) => snvs
                .iter()
                .map(|&snv| {
                    if snv < read.start_offset || snv > read.end_offset {
                        None
                    } else {
                        Some(mm_offsets.binary_search(&snv).is_ok())
                    }
                })
                .collect(),
            None => vec![None; snvs.len()],
        })
        .collect();

    profiles
}

fn call_snvs(region: (i32, i32), reads: &[HiFiRead], min_freq: f64) -> Vec<i32> {
    let offset_counts = reads
        .iter()
        .filter_map(|r| r.mismatch_offsets.as_ref())
        .flatten()
        .filter(|offset| region.0 <= **offset && **offset <= region.1)
        .counts();

    let total_reads = reads.len() as f64;
    offset_counts
        .into_iter()
        .filter(|(_, count)| *count as f64 / total_reads >= min_freq)
        .map(|(offset, _)| *offset)
        .sorted()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    fn make_read(encoding: &str) -> HiFiRead {
        let seq_start = encoding
            .find(|c: char| ['A', 'T', 'G', 'C'].contains(&c))
            .unwrap();
        let seq_end = encoding
            .rfind(|c: char| ['A', 'T', 'G', 'C'].contains(&c))
            .unwrap()
            + 1;
        let bases = encoding.as_bytes()[seq_start..seq_end].to_vec();
        let quals = "(".repeat(bases.len()).as_bytes().to_vec();
        let mismatches = encoding
            .as_bytes()
            .iter()
            .enumerate()
            .filter(|(_, char)| **char == b'X')
            .map(|(index, _)| {
                if index < seq_start {
                    index as i32 - seq_start as i32
                } else {
                    index as i32 - seq_end as i32
                }
            })
            .collect_vec();
        let start_offset = -(seq_start as i32);
        let end_offset = encoding.len() as i32 - seq_end as i32;

        HiFiRead {
            id: "read".to_string(),
            is_reverse: false,
            bases,
            quals,
            meth: None,
            read_qual: None,
            mismatch_offsets: Some(mismatches),
            start_offset,
            end_offset,
            cigar: None,
            hp_tag: None,
            mapq: 60,
        }
    }

    fn make_reads(encodings: &[&str]) -> Vec<HiFiRead> {
        encodings.iter().map(|e| make_read(e)).collect_vec()
    }

    #[test]
    fn if_het_snvs_then_genotype() {
        let reads = make_reads(&[
            "XX====TATATATA===X===",
            "XX=X==TATATATA===X===",
            "XX====TATATATATA=X=X===",
            "XX====TATATATATA=X=X===",
            "XX====TATATATATA=X=",
            "=TATATATA===X===",
        ]);
        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let result = genotype(&reads, &tr_seqs);

        let gt = Gt::from([
            TrSize {
                size: 8,
                ci: (8, 8),
            },
            TrSize {
                size: 10,
                ci: (10, 10),
            },
        ]);

        let alleles = vec!["TATATATA".to_string(), "TATATATATA".to_string()];
        let assignment = vec![0, 0, 1, 1, 1, 0];

        assert_eq!(result, Some((gt, alleles, assignment)));
    }

    #[test]
    fn if_hom_snvs_then_none() {
        let reads = make_reads(&[
            "XX====TATATATATA=X=X===",
            "XX====TATATATATA=X=X===",
            "XX====TATATATATA=X=X===",
            "XX====TATATATATA=X=X===",
        ]);

        let tr_seqs = reads
            .iter()
            .map(|r| std::str::from_utf8(&r.bases).unwrap())
            .collect_vec();
        let result = genotype(&reads, &tr_seqs);
        assert_eq!(result, None);
    }
}
