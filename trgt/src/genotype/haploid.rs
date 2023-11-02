use super::{Gt, TrSize};

pub fn genotype(sizes: &[usize], counts: &[usize]) -> Gt {
    let mut gts_and_penalties = Vec::new();
    for allele in sizes.iter() {
        let penalty = calc_gt_penalty(*allele, sizes, counts);
        gts_and_penalties.push((allele, penalty));
    }

    gts_and_penalties.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let size = gts_and_penalties.first().unwrap().0;
    let ci = (*sizes.iter().min().unwrap(), *sizes.iter().max().unwrap());

    let allele = TrSize::new(*size, ci);
    let mut gt = Gt::new();
    gt.push(allele);
    gt
}

fn calc_gt_penalty(allele: usize, sizes: &[usize], counts: &[usize]) -> f64 {
    let mut penalty = 0.0;

    for (size, count) in sizes.iter().zip(counts) {
        let term = if *size != allele {
            10.0 + 2.0 * allele.abs_diff(*size) as f64
        } else {
            0.0
        };
        penalty += term * (*count as f64);
    }

    penalty
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_tr() {
        let sizes = vec![3];
        let counts = vec![3];
        let gt = genotype(&sizes, &counts);

        let allele = TrSize::new(3, (3, 3));
        let mut expected = Gt::new();
        expected.push(allele);
        assert_eq!(gt, expected);
    }

    #[test]
    fn mosaic_tr() {
        let sizes = vec![10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
        let counts = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        let gt = genotype(&sizes, &counts);

        let allele = TrSize::new(50, (10, 100));
        let mut expected = Gt::new();
        expected.push(allele);
        assert_eq!(gt, expected);
    }

    #[test]
    fn tr_with_outliers() {
        let sizes = vec![10, 50];
        let counts = vec![4, 2];
        let gt = genotype(&sizes, &counts);

        let allele = TrSize::new(10, (10, 50));
        let mut expected = Gt::new();
        expected.push(allele);
        assert_eq!(gt, expected);
    }
}
