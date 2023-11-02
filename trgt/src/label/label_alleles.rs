use log::error;

use super::label_motif::label_motif;
use super::label_with_regexp::label_with_regexp;
use super::{label_with_hmm, struc::*};
use crate::label::Annotation;
use crate::locus::Locus;

pub fn label_alleles(locus: &Locus, alleles: &Vec<String>) -> Vec<Annotation> {
    let tr_type = get_tr_type(&locus.struc);

    match tr_type {
        TrType::RegExp => label_with_regexp(locus, alleles),
        TrType::SingleMotif => label_motif(&locus.struc, alleles),
        TrType::Hmm => label_with_hmm(locus, alleles),
    }
}

fn get_tr_type(encoding: &str) -> TrType {
    if encoding.matches('(').count() == 1 {
        TrType::SingleMotif
    } else if encoding.contains('(') {
        TrType::RegExp
    } else if encoding.contains('<') {
        TrType::Hmm
    } else {
        error!("Unexpected TR encoding {}", encoding);
        std::process::exit(1);
    }
}
