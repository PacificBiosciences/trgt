use super::read::Betas;

#[derive(Debug, Clone)]
pub struct AlleleAlign {
    pub seq: Align,
    pub reads: Vec<(Align, Betas)>,
}

pub type Align = Vec<AlignSeg>;

#[derive(Debug, Clone)]
pub struct AlignSeg {
    pub width: usize,
    pub op: AlignOp,
    pub seg_type: SegType,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AlignOp {
    Match,
    Subst,
    Ins,
    Del,
}

#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq)]
pub enum SegType {
    // The integer is the motif index or motifs.len() for unsegmented regions
    Tr(usize),
    LeftFlank,
    RightFlank,
}
