use crate::cluster::consensus;
use crate::cluster::math::median;
use crate::genotype::{Gt, TrSize};
use arrayvec::ArrayVec;
use bio::alignment::distance::levenshtein;
use bio::alignment::{pairwise::*, Alignment};
use itertools::Itertools;
use kodama::{linkage, Method};

pub fn genotype(seqs: &Vec<&[u8]>, trs: &[&str]) -> (Gt, Vec<String>, Vec<i32>) {
    const MAX_SEQS: usize = 50; // max number of sequences to do all-pairs Levenshtein
    let mut groups = cluster(seqs, MAX_SEQS);
    groups.sort_by_key(|a| a.len());

    let group1 = groups.pop().unwrap();
    let seqs1 = group1.iter().map(|i| trs[*i]).collect_vec();
    let backbone = *seqs1.first().unwrap();
    let aligns1 = align(backbone, &seqs1);
    let allele1 = consensus::repair_consensus(backbone, &seqs1, &aligns1);

    let size1 = TrSize {
        size: allele1.len(),
        ci: get_ci(&seqs1),
    };

    if groups.is_empty() {
        let gt = ArrayVec::from([size1.clone(), size1]);
        let alleles = vec![allele1.clone(), allele1];
        let classification = vec![0_i32; seqs.len()];
        return (gt, alleles, classification);
    }

    let group2 = groups.pop().unwrap();

    let tiny_group2 = group1.len() >= 10 && group2.len() == 1;
    let group1_len = median(&group1.iter().map(|i| trs[*i].len() as i32).collect_vec()).unwrap();
    let group2_len = median(&group2.iter().map(|i| trs[*i].len() as i32).collect_vec()).unwrap();
    let delta = (group1_len - group2_len).abs();
    let group1_frac = group1.len() as f64 / (group1.len() as f64 + group2.len() as f64);
    let small_group2 = delta <= 100.0 && group1_frac >= 0.80;

    if tiny_group2 || small_group2 {
        let gt = ArrayVec::from([size1.clone(), size1]);
        let alleles = vec![allele1.clone(), allele1];
        let classification = vec![0_i32; seqs.len()];
        return (gt, alleles, classification);
    }

    let seqs2 = group2.iter().map(|i| trs[*i]).collect_vec();
    let backbone = *seqs2.first().unwrap();
    let aligns2 = align(backbone, &seqs2);
    let allele2 = consensus::repair_consensus(backbone, &seqs2, &aligns2);

    let mut classifications = vec![2_i32; seqs.len()];

    // for reads already clustered, no need to do edit distance
    for seq_index in group1 {
        classifications[seq_index] = 0_i32;
    }
    for seq_index in group2 {
        classifications[seq_index] = 1_i32;
    }

    // assign reads that weren't used for clustering to the closest allele
    if seqs.len() > MAX_SEQS {
        let mut tie_breaker = 1;
        let a1 = &allele1.as_bytes();
        let a2 = &allele2.as_bytes();
        for (tr, classification) in trs.iter().zip(classifications.iter_mut()) {
            if *classification == 2 {
                //no allele assigned
                let tr = &tr.as_bytes();
                let dist1 = levenshtein(tr, a1);
                let dist2 = levenshtein(tr, a2);
                *classification = match dist1.cmp(&dist2) {
                    std::cmp::Ordering::Less => 0,
                    std::cmp::Ordering::Greater => 1,
                    std::cmp::Ordering::Equal => {
                        tie_breaker = (tie_breaker + 1) % 2;
                        tie_breaker
                    }
                };
            }
        }
    }
    let size2 = TrSize {
        size: allele2.len(),
        ci: get_ci(&seqs2),
    };
    let mut gt = Gt::new();
    gt.push(size1);
    gt.push(size2);
    (gt, vec![allele1, allele2], classifications)
}

pub fn cluster(seqs: &Vec<&[u8]>, max_seqs: usize) -> Vec<Vec<usize>> {
    if seqs.is_empty() {
        return Vec::new();
    }

    if seqs.len() == 1 {
        return vec![vec![0]];
    }

    if seqs.len() == 2 {
        return vec![vec![0], vec![1]];
    }

    let (seqs, orig_index) = if seqs.len() <= max_seqs {
        (seqs.clone(), (0..seqs.len()).collect::<Vec<usize>>())
    } else {
        // uniformly pick sequences by length distribution
        let num_reads = seqs.len();
        log::warn!(
            "Subsampling {} / {} reads for sequence-based clustering",
            max_seqs,
            num_reads
        );
        let mut ret: Vec<&[u8]> = Vec::new();
        let mut orig_index = vec![0; max_seqs];
        ret.reserve(max_seqs);
        let mut fast: f64 = 0.0;
        let step = (num_reads as f64) / (max_seqs as f64);

        for item in orig_index.iter_mut().take(max_seqs) {
            let ind = fast.floor() as usize;
            ret.push(seqs[ind]);
            *item = ind;
            fast += step;
        }

        (ret, orig_index)
    };

    let mut dists = get_dist_matrix(&seqs);
    let dendrogram = linkage(&mut dists, seqs.len(), Method::Ward);
    let cutoff = dendrogram
        .steps()
        .iter()
        .map(|s| s.dissimilarity)
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap()
        * 0.7;

    let mut num_groups = 0;
    let num_nodes = 2 * seqs.len() - 1;
    let mut membership = vec![None; num_nodes];

    for (cluster_index, step) in dendrogram.steps().iter().enumerate().rev() {
        let cluster = cluster_index + seqs.len();
        if step.dissimilarity <= cutoff {
            if membership[cluster].is_none() {
                membership[cluster] = Some(num_groups);
                num_groups += 1;
            }

            membership[step.cluster1] = membership[cluster];
            membership[step.cluster2] = membership[cluster];
        }
    }

    let mut groups = Vec::new();
    groups.reserve(seqs.len());

    for group in membership.into_iter().take(seqs.len()) {
        if let Some(group) = group {
            groups.push(group);
        } else {
            groups.push(num_groups);
            num_groups += 1;
        }
    }

    let mut seqs_by_group = vec![Vec::new(); num_groups];
    for (seq_index, group) in groups.iter().enumerate() {
        seqs_by_group[*group].push(orig_index[seq_index]);
    }

    seqs_by_group
}

fn get_ci(seqs: &[&str]) -> (usize, usize) {
    let min_val = seqs.iter().map(|x| x.len()).min().unwrap();
    let max_val = seqs.iter().map(|x| x.len()).max().unwrap();
    (min_val, max_val)
}

fn get_dist_matrix(seqs: &Vec<&[u8]>) -> Vec<f64> {
    let dist_len = seqs.len() * (seqs.len() - 1);
    let dist_len = (dist_len as f64 / 2.0) as usize;
    let mut dists = Vec::new();
    dists.reserve(dist_len);
    for index1 in 0..seqs.len() {
        for index2 in index1 + 1..seqs.len() {
            let dist = levenshtein(seqs[index1], seqs[index2]);
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
