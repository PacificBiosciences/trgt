use crate::align::AlignInfo;
use crate::locus::{Allele, BaseLabel, Locus};
use crate::pipe_plot::*;
use crate::struc::RegionLabel;

pub fn plot_genotype(
    locus: &Locus,
    params: &DisplayParams,
    genotype: &[Allele],
    aligns: &[AlignInfo],
) -> PipePlot {
    let tick_spacing = get_tick_spacing(locus, genotype);

    let mut allele_plots = Vec::new();
    for (allele_index, allele) in genotype.iter().enumerate() {
        let allele_aligns: Vec<&AlignInfo> = aligns
            .iter()
            .filter(|info| info.allele_index == allele_index)
            .collect();
        allele_plots.push(plot_allele(
            locus,
            allele,
            &allele_aligns,
            params.what_to_show == "meth",
            tick_spacing,
        ));
    }

    let mut labels = Vec::new();
    for motif in &locus.motifs {
        let color = get_color(locus, AlignOp::Match, &RegionLabel::Tr(0, 0, motif.clone()));
        labels.push((motif.clone(), color));
    }

    let height = 4;
    let legend = Legend { labels, height };
    PipePlot {
        panels: allele_plots,
        legend,
    }
}

fn plot_allele(
    locus: &Locus,
    allele: &Allele,
    aligns: &Vec<&AlignInfo>,
    show_betas: bool,
    tick_spacing: usize,
) -> PlotPanel {
    let labels = if show_betas {
        &allele.flank_labels
    } else {
        &allele.region_labels
    };

    let mut allele_plot = Vec::new();

    allele_plot.push(get_allele_pipe(
        locus,
        tick_spacing,
        &allele.region_labels,
        &allele.base_labels,
    ));
    for align in aligns {
        allele_plot.push(get_pipe(locus, labels, align, show_betas));
    }

    allele_plot
}

fn get_tick_spacing(_locus: &Locus, genotype: &[Allele]) -> usize {
    let max_motif_count = genotype
        .iter()
        .map(|a| {
            a.base_labels
                .iter()
                .filter(|l| **l == BaseLabel::MotifBound)
                .count()
        })
        .max()
        .unwrap();

    if max_motif_count <= 150 {
        1
    } else if max_motif_count <= 250 {
        5
    } else {
        10
    }
}
