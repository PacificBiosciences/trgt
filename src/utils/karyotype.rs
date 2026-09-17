use crate::utils::{Ploidy, Result};
use std::{collections::HashMap, fs, io::BufRead, path::Path};

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
                let path = Path::new(encoding);
                if !path.exists() {
                    if is_path_like(encoding) {
                        return Err(format!("Karyotype file does not exist: {encoding}"));
                    }
                    return Err(format!(
                        "Invalid karyotype: {encoding}. Expected XX, XY, or a path to an existing karyotype file"
                    ));
                }
                let file = fs::File::open(path)
                    .map_err(|e| format!("Karyotype file {}: {}", path.display(), e))?;
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
            let chrom = parts.next().ok_or(format!(
                "Missing chromosome/ploidy at line {}: {:?}",
                line_number + 1,
                line
            ))?;
            let ploidy_str = parts.next().ok_or(format!(
                "Missing chromosome/ploidy at line {}: {:?}",
                line_number + 1,
                line
            ))?;
            let ploidy = ploidy_str.parse::<Ploidy>().map_err(|e: String| {
                format!(
                    "Invalid ploidy at line {}: {:?}, {}",
                    line_number + 1,
                    line,
                    e
                )
            })?;

            if ploidies.contains_key(chrom) {
                Err(format!(
                    "Duplicate chromosome entry at line {}: {:?}",
                    line_number + 1,
                    line
                ))?
            } else {
                ploidies.insert(chrom.to_string(), ploidy);
            }
        }

        if ploidies.is_empty() {
            return Err("Karyotype file contains no chromosome ploidy entries".to_string());
        }

        Ok(Self {
            ploidy: PloidyInfo::Custom(ploidies),
        })
    }

    /// Returns how the karyotype was given, as a single line that two runs can be compared on.
    pub fn encoding(&self) -> String {
        match &self.ploidy {
            PloidyInfo::PresetXX => "XX".to_string(),
            PloidyInfo::PresetXY => "XY".to_string(),
            PloidyInfo::Custom(ploidies) => {
                let mut entries: Vec<String> = ploidies
                    .iter()
                    .map(|(chrom, ploidy)| format!("{chrom}:{ploidy:?}"))
                    .collect();
                entries.sort();
                entries.join(",")
            }
        }
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

fn is_path_like(encoding: &str) -> bool {
    encoding.contains(std::path::MAIN_SEPARATOR)
        || encoding.starts_with('.')
        || encoding.starts_with('~')
        || Path::new(encoding).extension().is_some()
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
        let error = Karyotype::from_reader(reader).unwrap_err();
        assert_eq!(error, "Missing chromosome/ploidy at line 2: \"chr2\"");
    }

    #[test]
    fn test_karyotype_from_reader_invalid_ploidy() {
        let data = "\
chr1 2\n\
chr2 3\n\
chrX 1\n";
        let reader = std::io::Cursor::new(data);
        let error = Karyotype::from_reader(reader).unwrap_err();
        assert_eq!(
            error,
            "Invalid ploidy at line 2: \"chr2 3\", must be set to 0, 1, or 2"
        );
    }

    #[test]
    fn test_karyotype_from_reader_duplicate_chromosome() {
        let data = "\
chr1 2\n\
chr1 1\n";
        let reader = std::io::Cursor::new(data);
        let error = Karyotype::from_reader(reader).unwrap_err();
        assert_eq!(error, "Duplicate chromosome entry at line 2: \"chr1 1\"");
    }

    #[test]
    fn test_karyotype_from_reader_empty_file() {
        let data = "";
        let reader = std::io::Cursor::new(data);
        let error = Karyotype::from_reader(reader).unwrap_err();
        assert_eq!(
            error,
            "Karyotype file contains no chromosome ploidy entries"
        );
    }

    #[test]
    fn test_karyotype_new_rejects_unknown_preset() {
        let error = Karyotype::new("XXY").unwrap_err();
        assert_eq!(
            error,
            "Invalid karyotype: XXY. Expected XX, XY, or a path to an existing karyotype file"
        );
    }

    #[test]
    fn test_karyotype_new_reports_missing_path_like_file() {
        let error = Karyotype::new("./missing.karyotype").unwrap_err();
        assert_eq!(error, "Karyotype file does not exist: ./missing.karyotype");
    }
}
