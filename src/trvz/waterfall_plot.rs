use super::color::get_meth_colors;
use super::color::Color;
use super::color::ColorMap;
use super::read::project_betas;
use super::read::Betas;
use super::{align::Align, locus::Locus, read::Read};
use crate::trvz::align::AlignOp;
use crate::trvz::align::AlignSeg;
use crate::trvz::align::SegType;
use crate::trvz::read::Beta;
use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment as BioAlign;
use bio::alignment::AlignmentOperation as BioOp;
use itertools::Itertools;
use pipeplot::Band;
use pipeplot::Legend;
use pipeplot::Pipe;
use pipeplot::PipePlot;
use pipeplot::Seg;
use pipeplot::Shape;

pub fn plot_waterfall(
    locus: &Locus,
    what_to_show: &str,
    reads: &[Read],
    colors: &ColorMap,
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

    plot(locus, what_to_show, &reads, colors)
}

fn align(locus: &Locus, longest_read: usize, read: &Read) -> (Align, Vec<Beta>) {
    let lf_ref = locus.left_flank.as_bytes();
    let rf_ref = locus.right_flank.as_bytes();
    let lf_read = read.seq[..lf_ref.len()].as_bytes();
    let rf_read = read.seq[read.seq.len() - locus.right_flank.len()..].as_bytes();
    let lf_bio_align = get_flank_align(lf_ref, lf_read);
    let mut align = convert(&lf_bio_align, SegType::LeftFlank);
    // Placeholder for TR alignment
    let tr = &read.seq[locus.left_flank.len()..read.seq.len() - locus.right_flank.len()];
    align.extend(label_motifs(&locus.motifs, tr));
    // Add deletion that lines up right flanks
    align.push(AlignSeg {
        width: longest_read - read.seq.len(),
        op: AlignOp::Del,
        seg_type: SegType::RightFlank,
    });
    let rf_bio_align = get_flank_align(rf_ref, rf_read);
    align.extend(convert(&rf_bio_align, SegType::RightFlank));

    let mut proj_betas = Vec::new();

    let lf_betas = &read
        .betas
        .iter()
        .filter(|beta| beta.pos < lf_read.len())
        .cloned()
        .collect_vec();
    proj_betas.extend(project_betas(&lf_bio_align, lf_betas));

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
        project_betas(&rf_bio_align, rf_betas)
            .iter()
            .map(|beta| Beta {
                pos: beta.pos + lf_read.len() + tr.len() + longest_read - read.seq.len(),
                value: beta.value,
            }),
    );

    // let mut betas = project_betas(&lf_bio_align, lf_betas);
    assert_eq!(
        read.betas.len(),
        lf_betas.len() + tr_betas.len() + rf_betas.len()
    );

    (align, proj_betas)
}

fn get_flank_align(ref_seq: &[u8], read_seq: &[u8]) -> BioAlign {
    let likely_max_len = ref_seq.len() + 50;
    let mut aligner = get_aligner(likely_max_len);
    aligner.global(read_seq, ref_seq)
}

/// Convert a rust-bio alignment into an internal alignments
fn convert(bio_align: &BioAlign, seg_type: SegType) -> Align {
    assert!(bio_align.xstart == 0);
    assert!(bio_align.ystart == 0);

    let mut align = Vec::new();
    let mut ref_pos = 0;
    let mut seq_pos = 0;

    for (bio_op, group) in &bio_align.operations.iter().group_by(|op| *op) {
        let run_len = group.count();

        let align_seg = match *bio_op {
            BioOp::Match => AlignSeg {
                width: run_len,
                op: AlignOp::Match,
                seg_type,
            },
            BioOp::Subst => AlignSeg {
                width: run_len,
                op: AlignOp::Subst,
                seg_type,
            },
            BioOp::Del => AlignSeg {
                width: run_len,
                op: AlignOp::Del,
                seg_type,
            },
            BioOp::Ins => AlignSeg {
                width: 0,
                op: AlignOp::Ins,
                seg_type,
            },
            _ => panic!("No logic to handle {bio_op:?}"),
        };

        align.push(align_seg);

        ref_pos += match *bio_op {
            BioOp::Match => run_len,
            BioOp::Subst => run_len,
            BioOp::Del => run_len,
            BioOp::Ins => 0,
            _ => panic!("Unhandled operation {:?}", *bio_op),
        };

        seq_pos += match *bio_op {
            BioOp::Match => run_len,
            BioOp::Subst => run_len,
            BioOp::Del => 0,
            BioOp::Ins => run_len,
            _ => panic!("Unhandled operation {:?}", *bio_op),
        };
    }

    assert_eq!(bio_align.x_aln_len(), seq_pos);
    assert_eq!(bio_align.y_aln_len(), ref_pos);

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

    PipePlot { pipes, legend }
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
    for ((op, seg_type), group) in &align.iter().group_by(|a| (a.op.clone(), a.seg_type)) {
        let width = group.into_iter().map(|seg| seg.width).sum();
        grouped_align.push(AlignSeg {
            width,
            op,
            seg_type,
        });
    }
    grouped_align
}

type ScoreFunc = fn(u8, u8) -> i32;

fn score(a: u8, b: u8) -> i32 {
    if a == b {
        1i32
    } else {
        -1i32
    }
}

fn get_aligner<'a>(likely_max_len: usize) -> Aligner<&'a ScoreFunc> {
    Aligner::with_capacity(
        likely_max_len,
        likely_max_len,
        -5,
        -1,
        &(score as ScoreFunc),
    )
}
