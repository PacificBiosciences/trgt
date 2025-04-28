use super::params::get_meth_colors;
use super::params::Color;
use super::params::ColorMap;
use super::params::PlotParams;
use super::read::project_betas;
use super::read::Betas;
use super::scale::get_scale;
use super::{align::Align, locus::Locus, read::Read};
use crate::trvz::align::AlignOp;
use crate::trvz::align::AlignSeg;
use crate::trvz::align::SegType;
use crate::trvz::align_consensus::align_motifs;
use crate::trvz::read::Beta;
use crate::wfaligner::AlignmentScope;
use crate::wfaligner::MemoryModel;
use crate::wfaligner::WFAligner;
use crate::wfaligner::WfaAlign;
use crate::wfaligner::WfaOp;
use itertools::Itertools;
use pipeplot::{Band, FontConfig, Legend, Pipe, PipePlot, Seg, Shape};

pub fn plot_waterfall(
    locus: &Locus,
    what_to_show: &str,
    reads: &[Read],
    params: &PlotParams,
) -> PipePlot {
    let reads = reads
        .iter()
        .sorted_by(|r1, r2| r1.seq.len().cmp(&r2.seq.len()))
        .collect_vec();

    let longest_read = reads.iter().map(|r| r.seq.len()).max().unwrap();
    let reads = reads
        .iter()
        .map(|r| align(locus, longest_read, r))
        .collect_vec();

    plot(locus, what_to_show, &reads, params)
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
    let motif_encoding = &locus
        .motifs
        .iter()
        .map(|m| m.as_bytes().to_vec())
        .collect_vec();
    let motif_aligns = align_motifs(motif_encoding, tr);
    align.extend(motif_aligns);
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
    params: &PlotParams,
) -> PipePlot {
    let xpos = 0;
    let mut ypos = 0;
    let mut pipes = Vec::new();
    let pipe = get_scale(xpos, ypos, params.pipe_height, &reads.last().unwrap().0);
    pipes.push(pipe);
    ypos += 4;
    for (align, betas) in reads {
        let (colors, betas) = if what_to_show == "meth" {
            (get_meth_colors(&locus.motifs), betas.clone())
        } else {
            (params.colors.clone(), Vec::new())
        };
        let pipe = get_pipe(xpos, ypos, params.pipe_height, align, &betas, &colors);
        pipes.push(pipe);
        ypos += params.pipe_height + params.pipe_pad;
    }

    let mut labels = Vec::new();
    if what_to_show == "motifs" {
        for (index, motif) in locus.motifs.iter().enumerate() {
            let color = params.colors.get(&SegType::Tr(index)).unwrap().to_string();
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
        height: 4,
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
            let shape = match align_seg.op {
                AlignOp::Del => Shape::HLine,
                AlignOp::Ins => Shape::VLine,
                AlignOp::Match | AlignOp::Subst => Shape::Rect,
            };
            let color = match align_seg.op {
                AlignOp::Match => colors.get(&align_seg.seg_type).unwrap(),
                AlignOp::Subst => &Color::Gray,
                AlignOp::Del => &Color::LightGray,
                _ => &Color::Black,
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
