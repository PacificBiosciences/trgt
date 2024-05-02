use crate::utils::{Ploidy, Result};
use std::{collections::HashMap, fs, io::BufRead};

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
    pub fn new(encoding: &str) -> Result<Self> {
        let ploidy = match encoding {
            "XX" => PloidyInfo::PresetXX,
            "XY" => PloidyInfo::PresetXY,
            _ => {
                let file =
                    fs::File::open(encoding).map_err(|e| format!("File {}: {}", encoding, e))?;
                let reader = std::io::BufReader::new(file);
                return Self::from_reader(reader);
            }
        };
        Ok(Self { ploidy })
    }

    #[cfg(test)]
    pub fn new_for_test(ploidies: HashMap<String, Ploidy>) -> Self {
        Self {
            ploidy: PloidyInfo::Custom(ploidies),
        }
    }

    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        let mut ploidies = HashMap::new();

        for (line_number, line) in reader.lines().enumerate() {
            let line =
                line.map_err(|e| format!("Error reading line {}: {}", line_number + 1, e))?;

            let mut parts = line.split_whitespace();
            let chrom = parts.next().ok_or(
                format!("Missing chromosome/ploidy at line {}", line_number + 1).to_string(),
            )?;
            let ploidy_str = parts.next().ok_or(
                format!("Missing chromosome/ploidy at line {}", line_number + 1).to_string(),
            )?;
            let ploidy = ploidy_str.parse::<Ploidy>().map_err(|e: String| {
                format!("Invalid ploidy at line {}, {}", line_number + 1, e)
            })?;

            if ploidies.contains_key(chrom) {
                Err(format!(
                    "Duplicate chromosome entry at line {}: {}",
                    line_number + 1,
                    chrom
                ))?
            } else {
                ploidies.insert(chrom.to_string(), ploidy);
            }
        }

        Ok(Self {
            ploidy: PloidyInfo::Custom(ploidies),
        })
    }

    pub fn get_ploidy(&self, chrom: &str) -> Result<Ploidy> {
        match &self.ploidy {
            PloidyInfo::PresetXX => match chrom {
                "Y" | "chrY" => Ok(Ploidy::Zero),
                _ => Ok(Ploidy::Two),
            },
            PloidyInfo::PresetXY => match chrom {
                "X" | "chrX" | "Y" | "chrY" => Ok(Ploidy::One),
                _ => Ok(Ploidy::Two),
            },
            PloidyInfo::Custom(ploidies) => {
                if let Some(&ploidy) = ploidies.get(chrom) {
                    Ok(ploidy)
                } else {
                    Err(format!(
                        "Ploidy was not specified for chromosome: {}",
                        chrom
                    ))
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::Ploidy;
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

    #[test]
    fn test_karyotype_from_reader() {
        let data = "\
chr1 2\n\
chr2 1\n\
chrX 1\n\
chrY 0\n";
        let reader = std::io::Cursor::new(data);
        let karyotype = Karyotype::from_reader(reader).unwrap();
        assert_eq!(karyotype.get_ploidy("chr1").unwrap(), Ploidy::Two);
        assert_eq!(karyotype.get_ploidy("chr2").unwrap(), Ploidy::One);
        assert_eq!(karyotype.get_ploidy("chrX").unwrap(), Ploidy::One);
        assert_eq!(karyotype.get_ploidy("chrY").unwrap(), Ploidy::Zero);
    }

    #[test]
    fn test_karyotype_from_reader_incomplete_line() {
        let data = "\
chr1 2\n\
chr2\n\
chrX 1\n";
        let reader = std::io::Cursor::new(data);
        let result = Karyotype::from_reader(reader);
        assert!(result.is_err());
    }

    #[test]
    fn test_karyotype_from_reader_invalid_ploidy() {
        let data = "\
chr1 2\n\
chr2 3\n\
chrX 1\n";
        let reader = std::io::Cursor::new(data);
        let result = Karyotype::from_reader(reader);
        assert!(result.is_err());
    }

    #[test]
    fn test_karyotype_from_reader_empty_file() {
        let data = "";
        let reader = std::io::Cursor::new(data);
        let result = Karyotype::from_reader(reader);
        assert!(result.is_ok());
        let karyotype = result.unwrap();
        assert!(karyotype.get_ploidy("chr1").is_err());
    }
}
