use crate::trvz::{
    align::AlignInfo,
    locus::Locus,
    pipe_plot::{
        get_color, get_pipe, AlignOp, Beta, Color, DisplayParams, Legend, Pipe, PipePlot, PipeSeg,
        Shape,
    },
    struc::RegionLabel,
};
use itertools::Itertools;

pub fn plot_waterfall(locus: &Locus, params: &DisplayParams, aligns: &[AlignInfo]) -> PipePlot {
    let flank_len = locus.left_flank.len();
    let mut aligns = aligns.iter().collect_vec();
    aligns.sort_by_key(|a| get_tr_len(flank_len, a));

    let max_tr_len = aligns
        .iter()
        .map(|a| get_tr_len(flank_len, a))
        .max()
        .unwrap();
    let mut waterfall = Vec::new();
    for align in aligns {
        let labels = if params.what_to_show == "motifs" {
            label_motifs(locus, align)
        } else {
            label_flanks(locus, align)
        };
        let pipe = get_pipe(locus, &labels, align, params.what_to_show == "meth");
        let pipe = add_deletion(max_tr_len, &labels, pipe);
        waterfall.push(pipe);
    }

    let mut labels = Vec::new();
    for motif in &locus.motifs {
        let color = get_color(locus, AlignOp::Match, &RegionLabel::Tr(0, 0, motif.clone()));
        labels.push((motif.clone(), color));
    }

    let legend = Legend { labels, height: 4 };
    PipePlot {
        panels: vec![waterfall],
        legend,
    }
}

fn label_motifs(locus: &Locus, align: &AlignInfo) -> Vec<RegionLabel> {
    let mut motifs = locus.motifs.clone();
    motifs.sort_by_key(|b| std::cmp::Reverse(b.len()));

    let mut labels = Vec::new();
    let flank_len = locus.left_flank.len();

    let tr_start = flank_len;
    let tr_end = align.align.y_aln_len() - flank_len;
    labels.push(RegionLabel::Flank(0, tr_start));

    let tr_seq = &align.seq[tr_start..tr_end];
    let mut index = 0;
    let mut motif_by_base = Vec::new();
    while index < tr_seq.len() {
        let mut motif_found = false;
        for (motif_index, motif) in motifs.iter().enumerate() {
            if index + motif.len() > tr_seq.len() {
                break;
            }

            if &tr_seq[index..index + motif.len()] == motif {
                let mut piece = vec![motif_index; motif.len()];
                motif_by_base.append(&mut piece);
                motif_found = true;
                index += motif.len();
                break;
            }
        }

        if !motif_found {
            motif_by_base.push(motifs.len());
            index += 1;
        }
    }

    let mut tr_pos = tr_start;
    let groups = motif_by_base.iter().group_by(|rec| *rec);
    for (motif_index, group) in &groups {
        let group = group.collect_vec();
        if *motif_index == motifs.len() {
            labels.push(RegionLabel::Other(tr_pos, tr_pos + group.len()));
        } else {
            let motif = motifs[*motif_index].clone();
            labels.push(RegionLabel::Tr(tr_pos, tr_pos + group.len(), motif));
        }

        tr_pos += group.len();
    }

    assert_eq!(tr_pos, tr_end);

    labels.push(RegionLabel::Flank(tr_end, align.align.y_aln_len()));
    labels
}

fn label_flanks(locus: &Locus, align: &AlignInfo) -> Vec<RegionLabel> {
    let flank_len = locus.left_flank.len();
    let mut labels = Vec::new();
    let tr_start = flank_len;
    let tr_end = align.align.y_aln_len() - flank_len;
    labels.push(RegionLabel::Flank(0, tr_start));
    labels.push(RegionLabel::Other(tr_start, tr_end));
    labels.push(RegionLabel::Flank(tr_end, align.align.y_aln_len()));
    labels
}

fn get_tr_len(flank_len: usize, align: &AlignInfo) -> usize {
    align.align.ylen - 2 * flank_len
}

fn add_deletion(max_tr_len: usize, labels: &[RegionLabel], pipe: Pipe) -> Pipe {
    let flank_len = match labels.first().unwrap() {
        RegionLabel::Flank(start, end) => end - start,
        RegionLabel::Other(start, end) => end - start,
        RegionLabel::Seq(start, end) => end - start,
        RegionLabel::Tr(start, end, _) => end - start,
    };

    let tr_end = *match labels.last().unwrap() {
        RegionLabel::Flank(start, _) => start,
        RegionLabel::Other(start, _) => start,
        RegionLabel::Seq(start, _) => start,
        RegionLabel::Tr(start, _, _) => start,
    };

    let tr_len = tr_end - flank_len;
    let deletion_len = max_tr_len - tr_len;
    if deletion_len == 0 {
        return pipe;
    }

    let mut segs = Vec::new();
    let mut ref_pos = 0;
    for seg in pipe.segs {
        if ref_pos == tr_end {
            segs.push(PipeSeg {
                width: deletion_len as u32,
                color: Color::Gray,
                shape: Shape::HLine,
            });
        }
        ref_pos += seg.width as usize;

        segs.push(seg);
    }

    let mut betas = Vec::new();
    for beta in pipe.betas {
        if beta.pos <= tr_end {
            betas.push(beta);
        } else {
            betas.push(Beta {
                value: beta.value,
                pos: beta.pos + deletion_len,
            });
        }
    }

    Pipe {
        segs,
        betas,
        height: pipe.height,
        outline: Vec::new(),
        scale: Vec::new(),
    }
}

// add tests

#[cfg(test)]
mod tests {
    #[test]
    fn test_parent_origin_matrix_1() {
        let mut motifs = ["AC", "GTG", "AAAAAA", "GGGGG"];
        // motifs.sort_by(|a, b| b.len().cmp(&a.len()));
        motifs.sort_by_key(|b| std::cmp::Reverse(b.len()));
    }
}
