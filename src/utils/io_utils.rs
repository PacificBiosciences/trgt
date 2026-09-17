use crate::utils::Result;
use rust_htslib::bcf;
use std::path::{Path, PathBuf};

/// Returns the path of an output file, formed by appending `output_suffix` to `output_prefix`.
pub fn output_path(output_prefix: &Path, output_suffix: &str) -> PathBuf {
    let mut path = output_prefix.to_path_buf().into_os_string();
    path.push(format!(".{output_suffix}"));
    PathBuf::from(path)
}

pub fn create_writer<T, F>(output_prefix: &Path, output_suffix: &str, f: F) -> Result<T>
where
    F: FnOnce(&str) -> Result<T>,
{
    f(output_path(output_prefix, output_suffix).to_str().unwrap())
}

#[derive(Debug, Clone)]
pub enum OutputType {
    Vcf {
        is_uncompressed: bool,
        level: Option<u8>,
    },
    Bcf {
        is_uncompressed: bool,
        level: Option<u8>,
    },
}

impl Default for OutputType {
    fn default() -> Self {
        OutputType::Vcf {
            is_uncompressed: false,
            level: None,
        }
    }
}

impl OutputType {
    pub fn infer_from_path(path: &Path) -> Self {
        let path_lower = path.as_os_str().to_string_lossy().to_lowercase();
        if path_lower.ends_with(".bcf.gz") {
            OutputType::Bcf {
                is_uncompressed: false,
                level: None,
            }
        } else if path_lower.ends_with(".vcf.gz") || path_lower.ends_with(".vcf.bgz") {
            OutputType::Vcf {
                is_uncompressed: false,
                level: None,
            }
        } else if path_lower.ends_with(".bcf") {
            OutputType::Bcf {
                is_uncompressed: true,
                level: None,
            }
        } else {
            OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }
        }
    }

    pub fn writer_options(&self) -> (bool, bcf::Format) {
        match self {
            OutputType::Vcf {
                is_uncompressed, ..
            } => (*is_uncompressed, bcf::Format::Vcf),
            OutputType::Bcf {
                is_uncompressed, ..
            } => (*is_uncompressed, bcf::Format::Bcf),
        }
    }

    pub fn output_suffix(&self) -> &'static str {
        match self {
            OutputType::Vcf {
                is_uncompressed, ..
            } => {
                if *is_uncompressed {
                    "vcf"
                } else {
                    "vcf.gz"
                }
            }
            OutputType::Bcf {
                is_uncompressed, ..
            } => {
                if *is_uncompressed {
                    "bcf"
                } else {
                    "bcf.gz"
                }
            }
        }
    }
}

impl std::fmt::Display for OutputType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OutputType::Vcf {
                is_uncompressed, ..
            } => {
                if *is_uncompressed {
                    write!(f, "uncompressed VCF")
                } else {
                    write!(f, "compressed VCF")
                }
            }
            OutputType::Bcf {
                is_uncompressed, ..
            } => {
                if *is_uncompressed {
                    write!(f, "uncompressed BCF")
                } else {
                    write!(f, "compressed BCF")
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf;

    #[test]
    fn default_output_type_is_compressed_vcf() {
        assert!(matches!(
            OutputType::default(),
            OutputType::Vcf {
                is_uncompressed: false,
                level: None
            }
        ));
    }

    #[test]
    fn infer_from_path_recognizes_known_extensions() {
        assert!(matches!(
            OutputType::infer_from_path(Path::new("out.vcf")),
            OutputType::Vcf {
                is_uncompressed: true,
                level: None
            }
        ));
        assert!(matches!(
            OutputType::infer_from_path(Path::new("out.VCF.GZ")),
            OutputType::Vcf {
                is_uncompressed: false,
                level: None
            }
        ));
        assert!(matches!(
            OutputType::infer_from_path(Path::new("out.vcf.bgz")),
            OutputType::Vcf {
                is_uncompressed: false,
                level: None
            }
        ));
        assert!(matches!(
            OutputType::infer_from_path(Path::new("out.bcf")),
            OutputType::Bcf {
                is_uncompressed: true,
                level: None
            }
        ));
        assert!(matches!(
            OutputType::infer_from_path(Path::new("out.bcf.gz")),
            OutputType::Bcf {
                is_uncompressed: false,
                level: None
            }
        ));

        assert!(matches!(
            OutputType::infer_from_path(Path::new("out.unknown")),
            OutputType::Vcf {
                is_uncompressed: true,
                level: None
            }
        ));
    }

    #[test]
    fn writer_options_return_expected_htslib_format() {
        let (is_uncompressed, format) = OutputType::Vcf {
            is_uncompressed: true,
            level: None,
        }
        .writer_options();
        assert!(is_uncompressed);
        assert_eq!(format, bcf::Format::Vcf);

        let (is_uncompressed, format) = OutputType::Bcf {
            is_uncompressed: false,
            level: None,
        }
        .writer_options();
        assert!(!is_uncompressed);
        assert_eq!(format, bcf::Format::Bcf);
    }

    #[test]
    fn output_suffix_matches_type_and_compression() {
        assert_eq!(OutputType::default().output_suffix(), "vcf.gz");
        assert_eq!(
            OutputType::Vcf {
                is_uncompressed: true,
                level: None
            }
            .output_suffix(),
            "vcf"
        );
        assert_eq!(
            OutputType::Bcf {
                is_uncompressed: true,
                level: None
            }
            .output_suffix(),
            "bcf"
        );
        assert_eq!(
            OutputType::Bcf {
                is_uncompressed: false,
                level: None
            }
            .output_suffix(),
            "bcf.gz"
        );
    }
}
