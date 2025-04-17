use crate::wfaligner::{WfaAlign, WfaOp};

#[derive(Clone)]
pub struct Read {
    pub read_name: String,
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
pub fn project_betas(bio_align: &WfaAlign, betas: &Betas) -> Betas {
    if betas.is_empty() {
        return Vec::new();
    }
    let mut ref_pos = 0;
    let mut seq_pos = 0;
    let mut beta_index = 0;

    let mut proj_betas = Vec::new();
    for op in &bio_align.operations {
        let at_pos = betas[beta_index].pos == seq_pos;
        let is_visible = *op == WfaOp::Match || *op == WfaOp::Subst;
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
            WfaOp::Match => 1,
            WfaOp::Subst => 1,
            WfaOp::Del => 0,
            WfaOp::Ins => 1,
        };

        ref_pos += match *op {
            WfaOp::Match => 1,
            WfaOp::Subst => 1,
            WfaOp::Del => 1,
            WfaOp::Ins => 0,
        };
    }

    proj_betas
}
