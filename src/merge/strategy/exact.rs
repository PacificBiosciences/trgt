use crate::utils::Result;
use rust_htslib::bcf::record::GenotypeAllele;
use std::collections::{HashMap, HashSet};

#[allow(clippy::type_complexity)]
pub fn merge_exact(
    vcf_gts: Vec<Vec<Vec<GenotypeAllele>>>,
    sample_alleles: Vec<Vec<&[u8]>>,
) -> Result<(Vec<Vec<Vec<GenotypeAllele>>>, Vec<&[u8]>)> {
    let mut ref_allele: Option<&[u8]> = None;
    let mut all_alleles: HashSet<&[u8]> = HashSet::new();

    for sample_allele in sample_alleles.iter() {
        if !sample_allele.is_empty() {
            if let Some(existing_ref_allele) = &ref_allele {
                if existing_ref_allele != &sample_allele[0] {
                    return Err(format_args!(
                        "Reference alleles do not match: '{}' and '{}'",
                        String::from_utf8_lossy(existing_ref_allele),
                        String::from_utf8_lossy(sample_allele[0])
                    )
                    .to_string());
                }
            } else {
                ref_allele = Some(sample_allele[0]);
            }
            for allele in &sample_allele[1..] {
                all_alleles.insert(allele);
            }
        }
    }
    let ref_allele = ref_allele.ok_or_else(|| "No reference allele found".to_string())?;

    let mut sorted_alleles: Vec<&[u8]> = all_alleles.into_iter().collect();
    sorted_alleles.sort_by(|a, b| match a.len().cmp(&b.len()) {
        std::cmp::Ordering::Equal => a.cmp(b),
        other => other,
    });
    sorted_alleles.insert(0, ref_allele);

    let allele_to_index: HashMap<&[u8], usize> = sorted_alleles
        .iter()
        .enumerate()
        .map(|(idx, &allele)| (allele, idx))
        .collect();

    let mut out_sample_gts = Vec::new();
    for (i, vcf_gt) in vcf_gts.iter().enumerate() {
        let mut out_gt = Vec::new();
        for sample_gt in vcf_gt {
            let mut s_gt = Vec::new();
            for gt in sample_gt {
                match gt {
                    GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing => {
                        s_gt.push(*gt)
                    }
                    GenotypeAllele::Phased(pos) | GenotypeAllele::Unphased(pos) => {
                        let pos_usize: usize = (*pos)
                            .try_into()
                            .map_err(|_| format!("Index out of range: {}", pos))?;
                        let pos_converted = sample_alleles[i]
                            .get(pos_usize)
                            .and_then(|allele| allele_to_index.get(allele))
                            .ok_or_else(|| {
                                format!(
                                    "Index out of bounds or allele not found in index: {:?}",
                                    &sample_alleles[i].get(pos_usize)
                                )
                            })?;
                        let new_gt = match gt {
                            GenotypeAllele::Phased(_) => {
                                GenotypeAllele::Phased(*pos_converted as i32)
                            }
                            GenotypeAllele::Unphased(_) => {
                                GenotypeAllele::Unphased(*pos_converted as i32)
                            }
                            _ => unreachable!(),
                        };
                        s_gt.push(new_gt);
                    }
                }
            }
            out_gt.push(s_gt);
        }
        out_sample_gts.push(out_gt);
    }
    Ok((out_sample_gts, sorted_alleles))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_exact() {
        let sample_gts = vec![
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
            vec![vec![
                GenotypeAllele::Unphased(0),
                GenotypeAllele::Unphased(0),
            ]],
            vec![vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing,
            ]],
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
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

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, sample_alleles).unwrap();

        assert_eq!(
            sorted_alleles[0],
            b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
        );

        assert_eq!(sorted_alleles.len(), 6);

        assert_eq!(
            out_gts[0][0],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)]
        );
        assert_eq!(
            out_gts[1][0],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(3)]
        );
        assert_eq!(
            out_gts[2][0],
            vec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
        );
        assert_eq!(
            out_gts[3][0],
            vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing
            ]
        );
        assert_eq!(
            out_gts[4][0],
            vec![GenotypeAllele::Unphased(4), GenotypeAllele::Unphased(5)]
        );
    }

    #[test]
    fn test_merge_exact_phasing() {
        let sample_gts = vec![
            vec![vec![
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
            ]],
            vec![vec![GenotypeAllele::Phased(1), GenotypeAllele::Phased(2)]],
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

        let (out_gts, sorted_alleles) = merge_exact(sample_gts, sample_alleles).unwrap();

        assert_eq!(
            sorted_alleles[0],
            b"CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA"
        );

        assert_eq!(sorted_alleles.len(), 5);

        assert_eq!(
            out_gts[0][0],
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)]
        );
        assert_eq!(
            out_gts[1][0],
            vec![GenotypeAllele::Phased(3), GenotypeAllele::Phased(4)]
        );
    }
}
