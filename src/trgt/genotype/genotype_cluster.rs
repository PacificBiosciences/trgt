use super::{consensus, Gt, TrSize};
use crate::{
    commands::genotype::THREAD_WFA_ED,
    utils::{align, Ploidy},
    wfa_aligner::WFAligner,
};
use arrayvec::ArrayVec;
use itertools::Itertools;
use kodama::{linkage, Method};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn central_read(num_seqs: usize, group: &[usize], dists: &[f64]) -> usize {
    let group_size = group.len();
    if group_size <= 2 {
        return group[0];
    }

    let mut dist_sums = vec![0_f64; group_size];
    for i in 0..(group_size - 1) {
        for j in (i + 1)..group_size {
            let index1 = group[i];
            let index2 = group[j];
            /* dist_sums has the condensed distance matrix, element 0 is
             * dist(seq[0], seq[1]), element 1 is dist(seq[0], seq[2]), etc.
             * This is the "inverse" that finds the matrix index based on the
             * index of the two seqs being compared */
            let mat_index = num_seqs * index1 - index1 * (index1 + 3) / 2 + index2 - 1;
            dist_sums[i] += dists[mat_index];
            dist_sums[j] += dists[mat_index];
        }
    }

    dist_sums
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(index, _)| group[index])
        .unwrap()
}

pub fn make_consensus(
    num_seqs: usize,
    trs: &[&[u8]],
    dists: &[f64],
    group: &[usize],
) -> (String, TrSize) {
    let seqs = group
        .iter()
        .map(|&i| std::str::from_utf8(trs[i]).unwrap())
        .collect_vec();
    let backbone = std::str::from_utf8(trs[central_read(num_seqs, group, dists)]).unwrap();
    let aligns = align(backbone, &seqs);
    let allele = consensus::repair_consensus(backbone, &seqs, &aligns);
    let size = TrSize::new(allele.len(), get_ci(&seqs));
    (allele, size)
}

pub fn genotype(ploidy: Ploidy, trs: &[&str]) -> (Gt, Vec<String>, Vec<i32>) {
    let trs: Vec<&[u8]> = trs.iter().map(|x| x.as_bytes()).collect();
    let mut dists = get_dist_matrix(&trs);
    let num_seqs = trs.len();
    if ploidy == Ploidy::One || num_seqs == 1 {
        let group: Vec<usize> = (0..num_seqs).collect();
        let (allele, size) = make_consensus(num_seqs, &trs, &dists, &group);
        let classifications = vec![0; num_seqs];
        if ploidy == Ploidy::One {
            let gt = Gt::from(size);
            return (gt, vec![allele], classifications);
        }

        // one read, two alleles
        let gt = Gt::from([size.clone(), size]);
        return (gt, vec![allele.clone(), allele], classifications);
    }
    let mut groups = cluster(num_seqs, &mut dists);

    assert!(groups.len() >= 2);
    groups.sort_by_key(|a| a.len());

    let group1 = groups.pop().unwrap();
    let group2 = groups.pop().unwrap();

    let (allele1, size1) = make_consensus(num_seqs, &trs, &dists, &group1);
    let (allele2, size2) = make_consensus(num_seqs, &trs, &dists, &group2);

    // GS: this should be handled better, but for now we just check for small
    // differences in size and large differences in cov
    fn small_group_is_outlier(len1: usize, len2: usize, cov1: usize, cov2: usize) -> bool {
        const MIN_LEN_DIFF: usize = 100;
        const MIN_COV_RATIO: usize = 4;

        let min_cov = std::cmp::min(cov1, cov2);
        let max_cov = std::cmp::max(cov1, cov2);
        len1.abs_diff(len2) < MIN_LEN_DIFF && min_cov * MIN_COV_RATIO < max_cov
    }
    if small_group_is_outlier(allele1.len(), allele2.len(), group1.len(), group2.len()) {
        // redo the homozygous case
        let group1: Vec<usize> = (0..num_seqs).step_by(2).collect();
        let group2: Vec<usize> = (1..num_seqs).step_by(2).collect();
        let (allele1, size1) = make_consensus(num_seqs, &trs, &dists, &group1);
        let (allele2, size2) = make_consensus(num_seqs, &trs, &dists, &group2);
        let mut classifications: Vec<i32> = (0..num_seqs).map(|x| (x % 2) as i32).collect();
        let (gt, alleles) = if allele1.len() > allele2.len() {
            classifications = classifications.iter().map(|x| 1 - x).collect();
            (ArrayVec::from([size2, size1]), vec![allele2, allele1])
        } else {
            (ArrayVec::from([size1, size2]), vec![allele1, allele2])
        };

        return (gt, alleles, classifications);
    }

    let mut classifications = vec![2; num_seqs];
    for seq_index in group1 {
        classifications[seq_index] = 0;
    }
    for seq_index in group2 {
        classifications[seq_index] = 1;
    }

    // assign outlier reads (discarded in cluster()) to the closest consensus
    // avoid constant conversion
    let a1 = allele1.as_bytes();
    let a2 = allele2.as_bytes();

    // NOTE: Should this run in parallel?
    THREAD_WFA_ED.with(|aligner_cell| {
        let mut aligner = aligner_cell.borrow_mut();
        for i in 0..num_seqs {
            let mut tie_breaker = 1;
            if classifications[i] == 2 {
                let dist1 = get_dist(trs[i], a1, &mut aligner);
                let dist2 = get_dist(trs[i], a2, &mut aligner);
                if dist1 < dist2 {
                    classifications[i] = 0;
                } else if dist2 < dist1 {
                    classifications[i] = 1;
                } else {
                    tie_breaker = (tie_breaker + 1) % 2;
                    classifications[i] = tie_breaker;
                }
            }
        }
    });

    let (gt, alleles) = if allele1.len() > allele2.len() {
        classifications = classifications.iter().map(|x| 1 - x).collect();
        (Gt::from([size2, size1]), vec![allele2, allele1])
    } else {
        (Gt::from([size1, size2]), vec![allele1, allele2])
    };

    (gt, alleles, classifications)
}

pub fn cluster(num_seqs: usize, dists: &mut [f64]) -> Vec<Vec<usize>> {
    assert!(num_seqs >= 2);
    assert_eq!(num_seqs * (num_seqs - 1) / 2, dists.len());
    if num_seqs == 2 {
        return vec![vec![0], vec![1]];
    }

    let dendrogram = linkage(dists, num_seqs, Method::Ward);

    let steps = dendrogram.steps();
    let mut cutoff = 0.0;

    const MIN_SMALLER_FRAC: f64 = 0.01;

    // for high-coverage: at least 10 reads must be on the smaller allele
    const MIN_CLUSTER_SIZE: usize = 2;

    // choose whichever is most liberal
    let min_cluster_size = std::cmp::max(
        MIN_CLUSTER_SIZE,
        (MIN_SMALLER_FRAC * (num_seqs as f64)).round() as usize,
    );
    for step in steps.iter().rev() {
        let size1 = dendrogram.cluster_size(step.cluster1);
        let size2 = dendrogram.cluster_size(step.cluster2);
        let min_size = std::cmp::min(size1, size2);
        if min_size >= min_cluster_size {
            cutoff = step.dissimilarity - 0.0001;
            break;
        }
    }

    // homozygous: split reads across alleles equally
    if cutoff == 0.0 {
        return vec![
            (0..num_seqs).step_by(2).collect(),
            (1..num_seqs).step_by(2).collect(),
        ];
    }

    let mut num_groups = 0;
    let num_nodes = 2 * num_seqs - 1;
    let mut membership = vec![None; num_nodes];

    for (cluster_index, step) in steps.iter().enumerate().rev() {
        let cluster = cluster_index + num_seqs;
        if step.dissimilarity <= cutoff {
            if membership[cluster].is_none() {
                membership[cluster] = Some(num_groups);
                num_groups += 1;
            }

            membership[step.cluster1] = membership[cluster];
            membership[step.cluster2] = membership[cluster];
        }
    }

    let mut groups = Vec::with_capacity(num_seqs);
    for group in membership.into_iter().take(num_seqs) {
        if let Some(group) = group {
            groups.push(group);
        } else {
            groups.push(num_groups);
            num_groups += 1;
        }
    }

    let mut seqs_by_group = vec![Vec::new(); num_groups];
    for (seq_index, group) in groups.iter().enumerate() {
        seqs_by_group[*group].push(seq_index);
    }

    seqs_by_group
}

fn get_ci(seqs: &[&str]) -> (usize, usize) {
    let min_val = seqs.iter().map(|x| x.len()).min().unwrap();
    let max_val = seqs.iter().map(|x| x.len()).max().unwrap();
    (min_val, max_val)
}

// We skip ED in cases we already know it will be too costly to do so
const MAX_OPS: usize = 10000;
#[inline(always)]
fn get_dist(seq1: &[u8], seq2: &[u8], aligner: &mut WFAligner) -> f64 {
    let seq_diff = seq1.len().abs_diff(seq2.len()) as i32; // lower bound on ED
    let num_operations = seq1.len() * seq2.len();
    let dist = if num_operations > MAX_OPS {
        seq_diff
    } else {
        let _status = aligner.align_end_to_end(seq1, seq2);
        aligner.score()
    };
    (dist as f64).sqrt()
}

fn get_dist_matrix(trs: &[&[u8]]) -> Vec<f64> {
    let n = trs.len();
    if n < 2 {
        return Vec::new();
    }
    let dist_len = n * (n - 1) / 2;

    // Parallelize the outer loop and then reuse the aligner in the inner loop
    let results_per_thread: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            // Calculate the base index part for the current i once before the inner loop, index for pair (i, j) in condensed upper triangle matrix
            // This represents the number of elements before the start of row 'i' in the condensed matrix, adjusted so that adding 'j' gives the final index.
            let base_index_for_i = i * n - (i * (i + 1)) / 2 - (i + 1);
            let mut thread_results = Vec::new();
            THREAD_WFA_ED.with(|aligner_cell| {
                // Access thread local and borrow ONCE per thread's chunk of i values, the inner loop reuses the same aligner instance
                let mut aligner = aligner_cell.borrow_mut();
                for j in (i + 1)..n {
                    let dist = get_dist(trs[i], trs[j], &mut aligner);
                    let index = base_index_for_i + j;
                    thread_results.push((index, dist));
                }
            });
            thread_results
        })
        .collect();

    // Flatten the results, bound check commented out.
    let mut dists = vec![0.0; dist_len];
    for thread_batch in results_per_thread {
        for (index, dist) in thread_batch {
            dists[index] = dist;
            // if index < dist_len {
            //     dists[index] = dist;
            // } else {
            //     // This should never happen if the index calculation is correct
            //     log::error!("Calculated index {} out of bounds ({})", index, dist_len,);
            // }
        }
    }
    assert_eq!(dists.len(), dist_len);
    dists
}
