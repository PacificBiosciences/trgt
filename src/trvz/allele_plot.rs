use super::align::SegType;
use super::align_allele::get_allele_align;
use super::params::{get_meth_colors, Color, ColorMap, PlotParams};
use super::read::{Betas, Read};
use super::scale::get_scale;
use crate::trvz::align::{Align, AlignOp};
use crate::trvz::locus::Locus;
use itertools::Itertools;
use pipeplot::{Band, FontConfig, Legend, Pipe, PipePlot, Seg, Shape};

pub fn plot_alleles(
    locus: &Locus,
    what_to_show: &str,
    allele_seqs: &[String],
    reads: &[Read],
    params: PlotParams,
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

    let allele_height = 4;
    let xpos = 0;
    let mut ypos = 0;
    let mut pipes = Vec::new();
    for (allele_index, allele) in aligns_by_allele.iter().enumerate() {
        let pipe = get_scale(xpos, ypos, allele_height, &allele.seq);
        pipes.push(pipe);
        ypos += allele_height;
        let pipe = get_pipe(
            xpos,
            ypos,
            allele_height,
            &allele.seq,
            &Vec::new(),
            &params.colors,
            true,
        );
        pipes.push(pipe);
        ypos += allele_height + params.pipe_pad;

        // Add extra padding if pipes are bookended
        if params.pipe_pad == 0 {
            ypos += 1;
        }

        // TODO: Confirm that this allele / index correspondence is always correct
        for (align, betas) in &allele.reads {
            let (colors, betas) = if what_to_show == "meth" {
                (get_meth_colors(&locus.motifs), betas.clone())
            } else {
                (params.colors.clone(), Vec::new())
            };

            let pipe = get_pipe(
                xpos,
                ypos,
                params.pipe_height,
                align,
                &betas,
                &colors,
                false,
            );
            pipes.push(pipe);
            ypos += params.pipe_height + params.pipe_pad;
        }

        if allele_index + 1 != aligns_by_allele.len() {
            ypos += 7;
        }
    }

    let mut labels = Vec::new();
    for (index, motif) in locus.motifs.iter().enumerate() {
        let color = params.colors.get(&SegType::Tr(index)).unwrap().to_string();
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
        height: allele_height,
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
