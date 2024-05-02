use std::str::FromStr;

#[derive(Debug, Clone, Copy)]
pub enum Genotyper {
    Size,
    Cluster,
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
