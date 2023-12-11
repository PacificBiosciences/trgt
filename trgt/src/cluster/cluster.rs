use crate::cluster::consensus;
use crate::cluster::math::median;
use crate::genotype::{Gt, TrSize};
use arrayvec::ArrayVec;
use bio::alignment::distance::simd::bounded_levenshtein;
use bio::alignment::{pairwise::*, Alignment};
use itertools::Itertools;
use kodama::{linkage, Method};

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
    trs: &[&str],
    dists: &[f64],
    group: &[usize],
) -> (String, TrSize) {
    let seqs = group.iter().map(|&i| trs[i]).collect_vec();
    let backbone = trs[central_read(num_seqs, group, dists)];
    let aligns = align(backbone, &seqs);
    let allele = consensus::repair_consensus(backbone, &seqs, &aligns);
    let size = TrSize::new(allele.len(), get_ci(&seqs));

    (allele, size)
}

pub fn genotype(seqs: &Vec<&[u8]>, trs: &[&str]) -> (Gt, Vec<String>, Vec<i32>) {
    let mut dists = get_dist_matrix(seqs);
    let mut groups = cluster(seqs.len(), &mut dists);
    groups.sort_by_key(|a| a.len());

    let group1 = groups.pop().unwrap();
    let (allele1, size1) = make_consensus(seqs.len(), trs, &dists, &group1);

    let group2 = groups.pop().unwrap_or_default();

    const MAX_GROUP_RATIO: usize = 10;
    let is_homozygous = group2.is_empty() || group2.len() * MAX_GROUP_RATIO <= group1.len() || {
        // reject group2 if estimated consensus size difference is negligible
        // and group 1 has significantly more reads.
        let group1_len =
            median(&group1.iter().map(|&i| trs[i].len() as i32).collect_vec()).unwrap();
        let group2_len =
            median(&group2.iter().map(|&i| trs[i].len() as i32).collect_vec()).unwrap();
        let delta = (group1_len - group2_len).abs();
        let group1_frac = group1.len() as f64 / (group1.len() as f64 + group2.len() as f64);

        const MIN_GROUP_SIZE_DIFF: f32 = 100.0;
        const MAX_SIMILAR_GROUP_FRAC: f64 = 0.80;
        delta <= MIN_GROUP_SIZE_DIFF && group1_frac >= MAX_SIMILAR_GROUP_FRAC
    };

    if is_homozygous {
        let gt = ArrayVec::from([size1.clone(), size1]);
        let alleles = vec![allele1.clone(), allele1];

        // distribute reads across alleles "randomly"
        let classification = (0..seqs.len() as i32).map(|x| x % 2).collect::<Vec<i32>>();
        return (gt, alleles, classification);
    }

    let (allele2, size2) = make_consensus(seqs.len(), trs, &dists, &group2);
    let mut classifications = vec![2; seqs.len()];
    for seq_index in group1 {
        classifications[seq_index] = 0;
    }
    for seq_index in group2 {
        classifications[seq_index] = 1;
    }

    let (gt, alleles) = if allele1.len() > allele2.len() {
        classifications = classifications.iter().map(|x| 1 - x).collect();
        (ArrayVec::from([size2, size1]), vec![allele2, allele1])
    } else {
        (ArrayVec::from([size1, size2]), vec![allele1, allele2])
    };
    (gt, alleles, classifications)
}

pub fn cluster(num_seqs: usize, dists: &mut Vec<f64>) -> Vec<Vec<usize>> {
    if num_seqs == 0 {
        return Vec::new();
    }

    assert_eq!(num_seqs * (num_seqs - 1) / 2, dists.len());
    if num_seqs == 1 {
        return vec![vec![0]];
    }

    if num_seqs == 2 {
        return vec![vec![0], vec![1]];
    }

    let dendrogram = linkage(dists, num_seqs, Method::Ward);

    // last element is the last merge, which is the highest dissimilarity
    let cutoff = dendrogram.steps().last().unwrap().dissimilarity;
    let cutoff = cutoff * 0.7;

    let mut num_groups = 0;
    let num_nodes = 2 * num_seqs - 1;
    let mut membership = vec![None; num_nodes];

    for (cluster_index, step) in dendrogram.steps().iter().enumerate().rev() {
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

fn get_dist_matrix(seqs: &Vec<&[u8]>) -> Vec<f64> {
    let dist_len = seqs.len() * (seqs.len() - 1) / 2;
    let mut dists = Vec::with_capacity(dist_len);
    for (index1, seq1) in seqs.iter().enumerate() {
        for (_index2, seq2) in seqs.iter().enumerate().skip(index1 + 1) {
            let max_len = std::cmp::max(seq1.len(), seq2.len()) as u32;
            let min_len = std::cmp::min(seq1.len(), seq2.len()) as u32;
            let length_diff = max_len - min_len;

            // we'll skip ED in cases we already know it will be too costly to do so.
            const MAX_K: u32 = 500;
            let dist = if length_diff <= MAX_K {
                bounded_levenshtein(seq1, seq2, MAX_K).unwrap_or(MAX_K)
            } else {
                length_diff // lower bound on ED
            };
            dists.push(dist as f64);
        }
    }
    assert_eq!(dists.len(), dist_len);
    dists
}

fn align(backbone: &str, seqs: &[&str]) -> Vec<Alignment> {
    let mut aligner = Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });
    seqs.iter()
        .map(|seq| aligner.global(seq.as_bytes(), backbone.as_bytes()))
        .collect_vec()
}
