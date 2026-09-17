use std::str::FromStr;

#[derive(Debug, Clone, Copy)]
pub enum Genotyper {
    Size,
    Cluster,
}

impl std::fmt::Display for Genotyper {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Genotyper::Size => write!(f, "size"),
            Genotyper::Cluster => write!(f, "cluster"),
        }
    }
}

impl FromStr for Genotyper {
    type Err = &'static str;
    fn from_str(genotyper: &str) -> Result<Self, Self::Err> {
        match genotyper {
            "size" => Ok(Genotyper::Size),
            "cluster" => Ok(Genotyper::Cluster),
            _ => Err("Invalid genotyper"),
        }
    }
}
