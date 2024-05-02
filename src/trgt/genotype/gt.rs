use arrayvec::ArrayVec;

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

impl From<TrSize> for Gt {
    fn from(tr_size: TrSize) -> Self {
        let mut gt = Gt::new();
        gt.push(tr_size);
        gt
    }
}
