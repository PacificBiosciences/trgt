use super::{
    align::{Align, AlignOp, AlignSeg, SegType},
    color::{get_meth_colors, Color, ColorMap},
    locus::Locus,
    read::{project_betas, Beta, Betas, Read},
};
use crate::wfaligner::{AlignmentScope, MemoryModel, WFAligner, WfaAlign, WfaOp};
use itertools::Itertools;
use pipeplot::{Band, FontConfig, Legend, Pipe, PipePlot, Seg, Shape};

pub fn plot_waterfall(
    locus: &Locus,
    what_to_show: &str,
    reads: &[Read],
    colors: &ColorMap,
) -> PipePlot {
    let sorted_reads = reads.iter().sorted_by_key(|r| r.seq.len()).collect_vec();
    let longest_read = sorted_reads.last().map_or(0, |r| r.seq.len());
    let aligned_reads = sorted_reads
        .iter()
        .map(|r| align(locus, longest_read, r))
        .collect_vec();
    plot(locus, what_to_show, &aligned_reads, colors)
}

fn align(locus: &Locus, longest_read: usize, read: &Read) -> (Align, Vec<Beta>) {
    let lf_ref = locus.left_flank.as_bytes();
    let rf_ref = locus.right_flank.as_bytes();

    let lf_read = read.seq[..lf_ref.len()].as_bytes();
    let rf_read = read.seq[read.seq.len() - locus.right_flank.len()..].as_bytes();

    let lf_wfa_align = get_flank_align(lf_ref, lf_read);
    let mut align = convert(&lf_wfa_align, SegType::LeftFlank);
    // Placeholder for TR alignment
    let tr = &read.seq[locus.left_flank.len()..read.seq.len() - locus.right_flank.len()];
    align.extend(label_motifs(&locus.motifs, tr));
    // Add deletion that lines up right flanks
    let deletion_width = longest_read.saturating_sub(read.seq.len());
    // prevents adding a 0-width segment
    if deletion_width > 0 {
        align.push(AlignSeg {
            width: deletion_width,
            op: AlignOp::Del,
            seg_type: SegType::RightFlank,
        });
    }

    let rf_wfa_align = get_flank_align(rf_ref, rf_read);
    align.extend(convert(&rf_wfa_align, SegType::RightFlank));

    let mut proj_betas = Vec::new();

    let lf_betas = &read
        .betas
        .iter()
        .filter(|beta| beta.pos < lf_read.len())
        .cloned()
        .collect_vec();
    proj_betas.extend(project_betas(&lf_wfa_align, lf_betas));

    let tr_betas = &read
        .betas
        .iter()
        .filter(|beta| lf_read.len() <= beta.pos && beta.pos < lf_read.len() + tr.len())
        .map(|beta| Beta {
            pos: beta.pos - lf_read.len(),
            value: beta.value,
        })
        .collect_vec();
    proj_betas.extend(tr_betas.iter().map(|beta| Beta {
        value: beta.value,
        pos: beta.pos + lf_read.len(),
    }));

    let rf_betas = &read
        .betas
        .iter()
        .filter(|beta| lf_read.len() + tr.len() <= beta.pos)
        .map(|beta| Beta {
            pos: beta.pos - lf_read.len() - tr.len(),
            value: beta.value,
        })
        .collect_vec();

    proj_betas.extend(
        project_betas(&rf_wfa_align, rf_betas)
            .iter()
            .map(|beta| Beta {
                pos: beta.pos + lf_read.len() + tr.len() + longest_read - read.seq.len(),
                value: beta.value,
            }),
    );

    assert_eq!(
        read.betas.len(),
        lf_betas.len() + tr_betas.len() + rf_betas.len()
    );

    (align, proj_betas)
}

fn get_flank_align(ref_seq: &[u8], read_seq: &[u8]) -> WfaAlign {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(2, 5, 1)
        .build();
    let _status = aligner.align_end_to_end(ref_seq, read_seq);
    aligner.get_alignment()
}

/// Convert a rust-wfa alignment into an internal alignments
fn convert(wfa_align: &WfaAlign, seg_type: SegType) -> Align {
    assert!(wfa_align.xstart == 0, "WFA alignment xstart should be 0");
    assert!(wfa_align.ystart == 0, "WFA alignment ystart should be 0");

    let mut align = Vec::new();
    let mut ref_pos = 0;
    let mut seq_pos = 0;

    for (wfa_op, group) in &wfa_align.operations.iter().chunk_by(|op| *op) {
        let run_len = group.count();

        let align_seg = match *wfa_op {
            WfaOp::Match => AlignSeg {
                width: run_len,
                op: AlignOp::Match,
                seg_type,
            },
            WfaOp::Subst => AlignSeg {
                width: run_len,
                op: AlignOp::Subst,
                seg_type,
            },
            WfaOp::Del => AlignSeg {
                width: run_len,
                op: AlignOp::Del,
                seg_type,
            },
            WfaOp::Ins => AlignSeg {
                width: 0,
                op: AlignOp::Ins,
                seg_type,
            },
        };

        align.push(align_seg);

        ref_pos += match *wfa_op {
            WfaOp::Match | WfaOp::Subst | WfaOp::Del => run_len,
            WfaOp::Ins => 0,
        };

        seq_pos += match *wfa_op {
            WfaOp::Match | WfaOp::Subst | WfaOp::Ins => run_len,
            WfaOp::Del => 0,
        };
    }

    assert_eq!(
        wfa_align.ylen, seq_pos,
        "Sequence length mismatch after conversion"
    );
    assert_eq!(
        wfa_align.xlen, ref_pos,
        "Reference length mismatch after conversion"
    );

    align
}

pub fn plot(
    locus: &Locus,
    what_to_show: &str,
    reads: &[(Align, Vec<Beta>)],
    colors: &ColorMap,
) -> PipePlot {
    let height = 4;
    let xpos = 0;
    let mut ypos = 0;
    let mut pipes = Vec::new();
    for (align, betas) in reads {
        let (colors, betas) = if what_to_show == "meth" {
            (get_meth_colors(&locus.motifs), betas.clone())
        } else {
            (colors.clone(), Vec::new())
        };
        let pipe = get_pipe(xpos, ypos, height, align, &betas, &colors);
        pipes.push(pipe);
        ypos += 5;
    }

    let mut labels = Vec::new();
    if what_to_show == "motifs" {
        for (index, motif) in locus.motifs.iter().enumerate() {
            let color = colors.get(&SegType::Tr(index)).unwrap().to_string();
            labels.push((motif.clone(), color));
        }
    } else {
        labels = vec![
            ("Methylated".to_string(), Color::Grad(1.0).to_string()),
            ("Unmethylated".to_string(), Color::Grad(0.0).to_string()),
        ];
    }
    ypos += 1;

    let legend = Legend {
        xpos,
        ypos,
        height,
        labels,
    };

    PipePlot {
        pipes,
        legend,
        font: FontConfig::default(),
    }
}

fn get_pipe(
    xpos: u32,
    ypos: u32,
    height: u32,
    align: &Align,
    betas: &Betas,
    colors: &ColorMap,
) -> Pipe {
    let segs = align
        .iter()
        .map(|align_seg| {
            let (shape, color) = match align_seg.op {
                AlignOp::Match => (Shape::Rect, colors.get(&align_seg.seg_type).unwrap()),
                AlignOp::Subst => (Shape::Rect, &Color::Gray),
                AlignOp::Del => (Shape::HLine, &Color::Black),
                AlignOp::Ins => (Shape::VLine, &Color::Black),
            };
            Seg {
                width: align_seg.width as u32,
                color: color.to_string(),
                shape,
            }
        })
        .collect();

    let bands = betas
        .iter()
        .map(|beta| Band {
            pos: beta.pos as u32,
            width: 2,
            color: Color::Grad(beta.value).to_string(),
        })
        .collect_vec();

    Pipe {
        xpos,
        ypos,
        height,
        segs,
        bands,
        outline: false,
    }
}

fn label_motifs(motifs: &[String], seq: &str) -> Align {
    // Preserve motif order after sorting
    let mut motifs = motifs.iter().enumerate().collect_vec();
    motifs.sort_by_key(|(_c, m)| std::cmp::Reverse(m.len()));

    let mut align = Vec::new();
    let mut pos = 0;
    while pos != seq.len() {
        let mut motif_found = false;
        for (motif_index, motif) in &motifs {
            if pos + motif.len() > seq.len() {
                continue;
            }

            if &seq[pos..pos + motif.len()] == *motif {
                align.push(AlignSeg {
                    width: motif.len(),
                    op: AlignOp::Match,
                    seg_type: SegType::Tr(*motif_index),
                });

                motif_found = true;
                pos += motif.len();
                break;
            }
        }

        if !motif_found {
            align.push(AlignSeg {
                width: 1,
                op: AlignOp::Match,
                seg_type: SegType::Tr(motifs.len()),
            });
            pos += 1;
        }
    }

    group(&align)
}

// TODO: factor out as a generic alignment operation
fn group(align: &Align) -> Align {
    let mut grouped_align = Vec::new();
    for ((op, seg_type), group) in &align.iter().chunk_by(|a| (a.op.clone(), a.seg_type)) {
        let width = group.into_iter().map(|seg| seg.width).sum();
        grouped_align.push(AlignSeg {
            width,
            op,
            seg_type,
        });
    }
    grouped_align
}
