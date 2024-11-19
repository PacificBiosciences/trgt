use super::align::{MotifBound, SegType};
use super::align_allele::get_allele_align;
use super::color::{get_meth_colors, Color, ColorMap};
use super::read::{Betas, Read};
use crate::trvz::align::{Align, AlignOp};
use crate::trvz::locus::Locus;
use itertools::Itertools;
use pipeplot::{Band, Legend, Pipe, PipePlot, Seg, Shape};

pub fn plot_alleles(
    locus: &Locus,
    what_to_show: &str,
    allele_seqs: &[String],
    reads: &[Read],
    colors: ColorMap,
) -> PipePlot {
    let aligns_by_allele = allele_seqs
        .iter()
        .enumerate()
        .map(|(index, allele_seq)| {
            let allele_reads = reads
                .iter()
                .filter(|r| r.allele == index as i32)
                .collect_vec();
            get_allele_align(locus, allele_seq, &allele_reads)
        })
        .collect_vec();
    let bounds_by_allele = aligns_by_allele
        .iter()
        .map(|align| &align.seq.1)
        .collect_vec();
    let tick_spacing = get_tick_spacing(&bounds_by_allele);

    let height = 4;
    let xpos = 0;
    let mut ypos = 0;
    let mut pipes = Vec::new();
    for (allele_index, allele) in aligns_by_allele.iter().enumerate() {
        let pipe = get_scales(
            xpos,
            ypos,
            height / 3,
            tick_spacing,
            locus.motifs.len(),
            &allele.seq.1,
        );
        pipes.push(pipe);
        ypos += height / 3;
        let pipe = get_pipe(
            xpos,
            ypos,
            height,
            &allele.seq.0,
            &Vec::new(),
            &colors,
            true,
        );
        pipes.push(pipe);
        ypos += 5;

        // TODO: Confirm that this allele / index correspondence is always correct
        for (align, betas) in &allele.reads {
            let (colors, betas) = if what_to_show == "meth" {
                (get_meth_colors(&locus.motifs), betas.clone())
            } else {
                (colors.clone(), Vec::new())
            };

            let pipe = get_pipe(xpos, ypos, height, align, &betas, &colors, false);
            pipes.push(pipe);
            ypos += 5;
        }

        if allele_index + 1 != aligns_by_allele.len() {
            ypos += 7;
        }
    }

    let mut labels = Vec::new();
    for (index, motif) in locus.motifs.iter().enumerate() {
        let color = colors.get(&SegType::Tr(index)).unwrap().to_string();
        labels.push((motif.clone(), color));
    }
    if what_to_show == "meth" {
        labels.push(("Methylated".to_string(), Color::Grad(1.0).to_string()));
        labels.push(("Unmethylated".to_string(), Color::Grad(0.0).to_string()));
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

fn get_scales(
    mut xpos: u32,
    ypos: u32,
    height: u32,
    tick_spacing: usize,
    motif_count: usize,
    bounds: &[MotifBound],
) -> Pipe {
    xpos += bounds.first().unwrap().start as u32;
    let mut segs = Vec::new();

    for (motif_index, group) in &bounds.iter().group_by(|bound| bound.motif_index) {
        let mut tick_index = 0;
        for bound in group {
            if tick_index % tick_spacing == 0 {
                let tick_label = if motif_index < motif_count {
                    Some(tick_index as u32)
                } else {
                    None
                };
                segs.push(Seg {
                    width: 0,
                    color: Color::Black.to_string(),
                    shape: Shape::Tick(tick_label),
                });
            }
            segs.push(Seg {
                width: (bound.end - bound.start) as u32,
                color: Color::Black.to_string(),
                shape: Shape::None,
            });
            tick_index += 1;
        }
        if tick_index % tick_spacing == 0 {
            let tick_label = if motif_index < motif_count {
                Some(tick_index as u32)
            } else {
                None
            };
            segs.push(Seg {
                width: 0,
                color: Color::Black.to_string(),
                shape: Shape::Tick(tick_label),
            });
        }
    }

    let mut culled_segs = Vec::new();
    for (is_tick, group) in &segs
        .into_iter()
        .group_by(|seg| matches!(seg.shape, Shape::Tick(_)))
    {
        if is_tick {
            culled_segs.push(group.last().unwrap());
        } else {
            culled_segs.extend(group);
        }
    }

    Pipe {
        xpos,
        ypos,
        height,
        segs: culled_segs,
        bands: Vec::new(),
        outline: false,
    }
}

fn get_pipe(
    xpos: u32,
    ypos: u32,
    height: u32,
    align: &Align,
    betas: &Betas,
    colors: &ColorMap,
    outline: bool,
) -> Pipe {
    let mut segs = Vec::new();
    for align_seg in align {
        let shape = match align_seg.op {
            AlignOp::Del => Shape::HLine,
            AlignOp::Ins => Shape::VLine,
            AlignOp::Match | AlignOp::Subst => Shape::Rect,
        };
        let color = if align_seg.op == AlignOp::Match {
            colors.get(&align_seg.seg_type).unwrap()
        } else if align_seg.op == AlignOp::Subst {
            &Color::Gray
        } else {
            &Color::Black
        };
        segs.push(Seg {
            width: align_seg.width as u32,
            color: color.to_string(),
            shape,
        });
    }

    let mut bands = Vec::new();

    for beta in betas {
        // Band width is 2 to cover the entire CpG
        let color = Color::Grad(beta.value);
        bands.push(Band {
            pos: beta.pos as u32,
            width: 2,
            color: color.to_string(),
        });
    }

    Pipe {
        xpos,
        ypos,
        height,
        segs,
        bands,
        outline,
    }
}

fn get_tick_spacing(bounds_by_allele: &[&Vec<MotifBound>]) -> usize {
    let max_motif_count = bounds_by_allele
        .iter()
        .map(|bounds| bounds.len())
        .max()
        .unwrap_or(1);

    match max_motif_count {
        0..=20 => 1,
        21..=100 => 5,
        101..=200 => 10,
        201..=500 => 20,
        501..=1000 => 50,
        1001..=1500 => 100,
        1501..=2000 => 200,
        2001..=3000 => 250,
        _ => 500,
    }
}
