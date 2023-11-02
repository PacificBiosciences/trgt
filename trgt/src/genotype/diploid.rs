use super::{Gt, TrSize};
use itertools::Itertools;
use std::cmp::{max, min};

pub fn genotype(sizes: &[usize], counts: &[usize]) -> Gt {
    let mut gts_and_penalties = Vec::new();
    for short_index in 0..sizes.len() {
        for long_index in short_index..sizes.len() {
            let gt = (sizes[short_index], sizes[long_index]);
            let penalty = calc_gt_penalty(&gt, sizes, counts);
            gts_and_penalties.push((gt, penalty));
        }
    }

    gts_and_penalties.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let top_gt = gts_and_penalties.first().unwrap().0;
    let mut short_size = std::cmp::min(top_gt.0, top_gt.1);
    let mut long_size = std::cmp::max(top_gt.0, top_gt.1);

    if short_size != long_size && sizes.len() >= 2 {
        let coverage = counts.iter().sum::<usize>();
        let length_hist = sizes
            .iter()
            .zip(counts.iter())
            .sorted_by(|a, b| b.1.cmp(a.1))
            .collect_vec();

        let top_frac = *length_hist[0].1 as f64 / coverage as f64;
        let max_len = *sizes.iter().max().unwrap();
        let min_len = *sizes.iter().min().unwrap();
        let range = max_len - min_len;

        // Reclassify the repeat as homozygous if the entire range of lengths
        // is <= 6bp and > 60% of reads support a repeat the same length
        if top_frac > 0.60 && range <= 6 {
            let top_len = *length_hist[0].0;
            short_size = top_len;
            long_size = top_len;
        }
    }

    let (short_ci, long_ci) = get_ci((short_size, long_size), sizes);

    let short_allele = TrSize {
        size: short_size,
        ci: short_ci,
    };

    let long_allele = TrSize {
        size: long_size,
        ci: long_ci,
    };

    Gt::from([short_allele, long_allele])
}

fn calc_gt_penalty(gt: &(usize, usize), sizes: &[usize], counts: &[usize]) -> f64 {
    let (short_allele, long_allele) = gt;

    let mut penalty = 0.0;

    // max_frac roughly corresponds to the fraction of discrepant reads that
    // don't agree with any allele; we decrease this parameter to 0.05 when
    // one allele is significantly expanded over the other to account for the
    // possibility of lower coverage in the longer allele
    let max_frac = if short_allele.abs_diff(*long_allele) <= 100 {
        0.25
    } else {
        0.05
    };

    for (size, count) in sizes.iter().zip(counts) {
        let short_term = if *size != *short_allele {
            10 + 2 * short_allele.abs_diff(*size)
        } else {
            0
        };

        let long_term = if *size != *long_allele {
            10 + 2 * long_allele.abs_diff(*size)
        } else {
            0
        };

        let term = min(short_term, long_term) as f64 + max_frac * max(short_term, long_term) as f64;
        penalty += term * (*count as f64);
    }

    penalty
}

fn get_ci(gt: (usize, usize), sizes: &[usize]) -> ((usize, usize), (usize, usize)) {
    let (short_size, long_size) = gt;
    let mut short_ci = (short_size, short_size);
    let mut long_ci = (long_size, long_size);

    for size in sizes {
        if size.abs_diff(short_size) <= size.abs_diff(long_size) {
            short_ci = (min(short_ci.0, *size), max(short_ci.1, *size));
        } else {
            long_ci = (min(long_ci.0, *size), max(long_ci.1, *size));
        }
    }

    let short_ci = (short_ci.0, short_ci.1);
    let long_ci = (long_ci.0, long_ci.1);

    (short_ci, long_ci)
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrayvec::ArrayVec;

    #[test]
    fn clean_het_tr() {
        let sizes = vec![3, 4];
        let counts = vec![3, 3];
        let gt = genotype(&sizes, &counts);

        let short_allele = TrSize::new(3, (3, 3));
        let long_allele = TrSize::new(4, (4, 4));
        assert_eq!(gt, ArrayVec::from([short_allele, long_allele]));
    }
}
