use bio::alignment::Alignment as BioAlign;
use bio::alignment::AlignmentOperation as BioOp;

#[derive(Clone)]
pub struct Read {
    pub seq: String,
    pub left_flank: usize,
    pub right_flank: usize,
    pub allele: i32,
    pub betas: Betas,
}

pub type Betas = Vec<Beta>;

#[derive(Debug, Clone)]
pub struct Beta {
    pub pos: usize,
    pub value: f64,
}

/// Convert betas from read coordinates to allele coordinates
/// according to the provided alignment between the read and
/// the allele consensus sequence
pub fn project_betas(bio_align: &BioAlign, betas: &Betas) -> Betas {
    if betas.is_empty() {
        return Vec::new();
    }
    let mut ref_pos = 0;
    let mut seq_pos = 0;
    let mut beta_index = 0;

    let mut proj_betas = Vec::new();
    for op in &bio_align.operations {
        let at_pos = betas[beta_index].pos == seq_pos;
        let is_visible = *op == BioOp::Match || *op == BioOp::Subst;
        if at_pos && is_visible {
            let beta = Beta {
                pos: ref_pos,
                value: betas[beta_index].value,
            };
            proj_betas.push(beta);
        }
        if at_pos {
            beta_index += 1;
        }
        if beta_index == betas.len() {
            break;
        }

        seq_pos += match *op {
            BioOp::Match => 1,
            BioOp::Subst => 1,
            BioOp::Del => 0,
            BioOp::Ins => 1,
            _ => panic!("Unhandled operation {:?}", *op),
        };

        ref_pos += match *op {
            BioOp::Match => 1,
            BioOp::Subst => 1,
            BioOp::Del => 1,
            BioOp::Ins => 0,
            _ => panic!("Unhandled operation {op:?}"),
        };
    }

    proj_betas
}
