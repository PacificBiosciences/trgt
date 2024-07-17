use rust_htslib::bcf::record::GenotypeAllele;
use std::collections::{HashMap, HashSet};

pub fn merge_exact(
    sample_gts: Vec<Vec<GenotypeAllele>>,
    sample_alleles: Vec<Vec<&[u8]>>,
) -> (Vec<Vec<GenotypeAllele>>, Vec<&[u8]>) {
    let mut ref_allele: Option<&[u8]> = None;
    let mut all_alleles: HashSet<&[u8]> = HashSet::new();

    for sample_allele in sample_alleles.iter() {
        if !sample_allele.is_empty() {
            if let Some(ref_allele) = &ref_allele {
                assert_eq!(
                    ref_allele, &sample_allele[0],
                    "Reference alleles do not match"
                );
            } else {
                ref_allele = Some(sample_allele[0]);
            }
            for allele in &sample_allele[1..] {
                all_alleles.insert(allele);
            }
        }
    }
    let ref_allele = ref_allele.expect("No reference allele found");

    let mut sorted_alleles: Vec<&[u8]> = all_alleles.into_iter().collect();
    sorted_alleles.sort_by_key(|a| a.len());
    sorted_alleles.insert(0, ref_allele);

    let allele_to_index: HashMap<&[u8], usize> = sorted_alleles
        .iter()
        .enumerate()
        .map(|(idx, &allele)| (allele, idx))
        .collect();

    let mut out_sample_gts: Vec<Vec<GenotypeAllele>> = Vec::new();
    for (i, sample_gt) in sample_gts.iter().enumerate() {
        let mut out_gt: Vec<GenotypeAllele> = Vec::new();
        for gt in sample_gt {
            match gt {
                GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing => out_gt.push(*gt),
                GenotypeAllele::Phased(pos) | GenotypeAllele::Unphased(pos) => {
                    let pos_usize: usize = (*pos).try_into().expect("Index out of range");
                    let pos_converted = allele_to_index[&sample_alleles[i][pos_usize]];
                    let new_gt = match gt {
                        GenotypeAllele::Phased(_) => GenotypeAllele::Phased(pos_converted as i32),
                        GenotypeAllele::Unphased(_) => {
                            GenotypeAllele::Unphased(pos_converted as i32)
                        }
                        _ => unreachable!(),
                    };
                    out_gt.push(new_gt);
                }
            }
        }
        out_sample_gts.push(out_gt);
    }
    (out_sample_gts, sorted_alleles)
}

#[allow(dead_code)]
fn vec_to_comma_separated_string(vec: Vec<&[u8]>) -> String {
    vec.into_iter()
        .map(|slice| String::from_utf8_lossy(slice).to_string())
        .collect::<Vec<String>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_exact() {
        let sample_gts = vec![
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)],
            vec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)],
            vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing,
            ],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)],
        ];

        let sample_alleles = vec![
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA".as_ref(),
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA".as_ref(),
            ],
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA".as_ref(),
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
            ],
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
            ],
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
            ],
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
                b"CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGCGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTG"
                    .as_ref(),
                b"CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGCGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTG"
                    .as_ref(),
            ],
        ];

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, sample_alleles);

        assert_eq!(
            sorted_alleles[0],
            b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
        );

        assert_eq!(sorted_alleles.len(), 6);

        assert_eq!(
            out_gts[0],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)]
        );
        assert_eq!(
            out_gts[1],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(3)]
        );
        assert_eq!(
            out_gts[2],
            vec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
        );
        assert_eq!(
            out_gts[3],
            vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing
            ]
        );
        assert_eq!(
            out_gts[4],
            vec![GenotypeAllele::Unphased(4), GenotypeAllele::Unphased(5)]
        );
    }

    #[test]
    fn test_merge_exact_phasing() {
        let sample_gts = vec![
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)],
            vec![GenotypeAllele::Phased(1), GenotypeAllele::Phased(2)],
        ];

        let sample_alleles = vec![
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA".as_ref(),
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA".as_ref(),
            ],
            vec![
                b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
                    .as_ref(),
                b"CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGCGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTG"
                    .as_ref(),
                b"CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGCGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTGGGGTG"
                    .as_ref(),
            ],
        ];

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, sample_alleles);

        assert_eq!(
            sorted_alleles[0],
            b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
        );

        assert_eq!(sorted_alleles.len(), 5);

        assert_eq!(
            out_gts[0],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)]
        );
        assert_eq!(
            out_gts[1],
            vec![GenotypeAllele::Phased(3), GenotypeAllele::Phased(4)]
        );
    }
}
