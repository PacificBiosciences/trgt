use crate::genotype::Ploidy;
use std::collections::HashMap;
use std::fs;

#[derive(Debug, PartialEq, Clone)]
pub struct Karyotype {
    ploidy: PloidyInfo,
}

#[derive(Debug, PartialEq, Clone)]
enum PloidyInfo {
    PresetXX,
    PresetXY,
    Custom(HashMap<String, Ploidy>),
}

impl Karyotype {
    pub fn new(encoding: &str) -> Result<Self, String> {
        let ploidy = match encoding {
            "XX" => PloidyInfo::PresetXX,
            "XY" => PloidyInfo::PresetXY,
            _ => return Self::from_file(encoding),
        };
        Ok(Self { ploidy })
    }

    #[cfg(test)]
    pub fn new_for_test(ploidies: HashMap<String, Ploidy>) -> Self {
        Self {
            ploidy: PloidyInfo::Custom(ploidies),
        }
    }

    fn from_file(path: &str) -> Result<Self, String> {
        let contents = fs::read_to_string(path).map_err(|e| format!("File {}: {}", path, e))?;

        let ploidies = contents
            .lines()
            .enumerate()
            .map(|(line_number, line)| {
                let mut parts = line.split_whitespace();
                let chrom = parts
                    .next()
                    .ok_or(format!("Missing chromosome at line {}", line_number).to_string())?;
                let ploidy_str = parts
                    .next()
                    .ok_or(format!("Missing ploidy at line {}", line_number).to_string())?;
                let ploidy = ploidy_str.parse().map_err(|e: String| {
                    format!("Invalid ploidy at line {}, {}", line_number, e)
                })?;
                Ok((chrom.to_string(), ploidy))
            })
            .collect::<Result<HashMap<_, _>, String>>()?;

        Ok(Self {
            ploidy: PloidyInfo::Custom(ploidies),
        })
    }

    pub fn get_ploidy(&self, chrom: &str) -> Result<Ploidy, String> {
        match &self.ploidy {
            PloidyInfo::PresetXX => match chrom {
                "Y" | "chrY" => Ok(Ploidy::Zero),
                _ => Ok(Ploidy::Two),
            },
            PloidyInfo::PresetXY => match chrom {
                "X" | "chrX" | "Y" | "chrY" => Ok(Ploidy::One),
                _ => Ok(Ploidy::Two),
            },
            PloidyInfo::Custom(ploidies) => ploidies
                .get(chrom)
                .copied()
                .ok_or_else(|| format!("Ploidy was not specified for chromosome: {}", chrom)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::Ploidy;
    use std::collections::HashMap;

    #[test]
    fn test_karyotype_custom() {
        let mut ploidies = HashMap::new();
        ploidies.insert("chr1".to_string(), Ploidy::Two);
        ploidies.insert("chr2".to_string(), Ploidy::One);

        let karyotype = Karyotype::new_for_test(ploidies);

        assert_eq!(karyotype.get_ploidy("chr1").unwrap(), Ploidy::Two);
        assert_eq!(karyotype.get_ploidy("chr2").unwrap(), Ploidy::One);
        assert!(karyotype.get_ploidy("chrX").is_err());
    }
}
