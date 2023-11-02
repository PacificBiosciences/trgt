use crate::{label::Annotation, reads::HiFiRead};
use arrayvec::ArrayVec;

#[derive(Debug)]
pub struct Allele {
    pub seq: String,
    pub annotation: Annotation,
    pub ci: (usize, usize),
    pub num_spanning: usize,
    pub meth: Option<f64>,
}

pub type Genotype = ArrayVec<Allele, 2>;

#[derive(Debug)]
pub struct LocusResult {
    pub genotype: Genotype,
    pub reads: Vec<HiFiRead>,
    pub tr_spans: Vec<(usize, usize)>,
    pub classification: Vec<i32>,
}

impl LocusResult {
    pub fn empty() -> LocusResult {
        LocusResult {
            genotype: Genotype::new(),
            reads: Vec::new(),
            tr_spans: Vec::new(),
            classification: Vec::new(),
        }
    }
}
