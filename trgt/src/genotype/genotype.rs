use super::consensus::repair_consensus;
use super::diploid;
use super::haploid;
use super::Gt;
use super::Ploidy;
use crate::locate::get_consensus;
use arrayvec::ArrayVec;
use itertools::Itertools;

pub fn genotype(ploidy: Ploidy, seqs: &Vec<&str>) -> (Gt, Vec<String>, Vec<i32>) {
    let (unique_lens, len_counts) = get_len_hist(seqs);

    let gt = match ploidy {
        Ploidy::Zero => panic!("Can't genotype repeats of zero ploidy"),
        Ploidy::One => haploid::genotype(&unique_lens, &len_counts),
        Ploidy::Two => diploid::genotype(&unique_lens, &len_counts),
    };

    let allele_lens: ArrayVec<usize, 2> = gt.iter().map(|a| a.size).collect();
    let (unique_seqs, counts) = get_seq_hist(seqs.clone());
    let mut alleles = get_consensus(allele_lens.clone(), &unique_seqs, &counts);
    let seqs_by_allele = split(allele_lens, &unique_seqs, &counts);

    let mut fixed_alleles = Vec::new();
    for (index, allele) in alleles.iter().enumerate() {
        let seqs = &seqs_by_allele[index].0;
        let counts = &seqs_by_allele[index].1;
        let fixed_allele = repair_consensus(allele, seqs, counts);
        fixed_alleles.push(fixed_allele);
    }
    alleles = fixed_alleles;

    let mut classifications = vec![0_i32; seqs.len()];
    let mut tie_breaker = 1;
    for (seq, classification) in seqs.iter().zip(classifications.iter_mut()) {
        if alleles.len() == 2 {
            let diff1 = seq.len().abs_diff(alleles[0].len());
            let diff2 = seq.len().abs_diff(alleles[1].len());
            *classification = match diff1.cmp(&diff2) {
                std::cmp::Ordering::Less => 0,
                std::cmp::Ordering::Greater => 1,
                std::cmp::Ordering::Equal => {
                    tie_breaker = (tie_breaker + 1) % 2;
                    tie_breaker
                }
            };
        }
    }
    // allele_seqs is expected to contain two elements
    if alleles.len() == 1 {
        alleles.push(alleles[0].clone());
    }

    (gt, alleles, classifications)
}

fn get_len_hist(seqs: &[&str]) -> (Vec<usize>, Vec<usize>) {
    let sorted_lens = seqs.iter().map(|seq| seq.len()).sorted().collect_vec();
    let mut unique_lens = Vec::new();
    let mut unique_len_counts = Vec::new();
    for (len, group) in &sorted_lens.iter().group_by(|l| **l) {
        unique_lens.push(len);
        unique_len_counts.push(group.count());
    }

    (unique_lens, unique_len_counts)
}

fn get_seq_hist(mut seqs: Vec<&str>) -> (Vec<&str>, Vec<usize>) {
    seqs.sort();
    let mut unique_seqs = Vec::new();
    let mut seq_counts = Vec::new();
    for (seq, seqs) in &seqs.iter().group_by(|s| **s) {
        let count = seqs.count();
        unique_seqs.push(seq);
        seq_counts.push(count);
    }

    (unique_seqs, seq_counts)
}

fn split<'a>(
    allele_lens: ArrayVec<usize, 2>,
    seqs: &'a [&'a str],
    counts: &'a [usize],
) -> Vec<(Vec<&'a str>, Vec<usize>)> {
    if allele_lens.len() == 1 {
        return vec![(seqs.to_vec(), counts.to_vec())];
    }
    let (al1, al2) = allele_lens.iter().collect_tuple().unwrap();
    let al1_seqs: Vec<&str> = seqs
        .iter()
        .filter(|s| s.len().abs_diff(*al1) <= s.len().abs_diff(*al2))
        .cloned()
        .collect_vec();

    let al1_counts = seqs
        .iter()
        .zip(counts)
        .filter(|rec| al1_seqs.contains(rec.0))
        .map(|rec| *rec.1)
        .collect_vec();

    let al2_seqs = seqs
        .iter()
        .filter(|s| s.len().abs_diff(*al2) < s.len().abs_diff(*al1))
        .cloned()
        .collect_vec();

    let al2_counts = seqs
        .iter()
        .zip(counts)
        .filter(|rec| al2_seqs.contains(rec.0))
        .map(|rec| *rec.1)
        .collect_vec();

    vec![(al1_seqs, al1_counts), (al2_seqs, al2_counts)]
}
