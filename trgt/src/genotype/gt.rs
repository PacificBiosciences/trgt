use arrayvec::ArrayVec;
use std::str::FromStr;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Ploidy {
    Zero,
    One,
    Two,
}

impl FromStr for Ploidy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "0" => Ok(Ploidy::Zero),
            "1" => Ok(Ploidy::One),
            "2" => Ok(Ploidy::Two),
            _ => Err("must be set to 0, 1, or 2".to_string()),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct TrSize {
    pub size: usize,
    pub ci: (usize, usize),
}

impl TrSize {
    pub fn new(size: usize, ci: (usize, usize)) -> TrSize {
        TrSize { size, ci }
    }
}

pub type Gt = ArrayVec<TrSize, 2>;
