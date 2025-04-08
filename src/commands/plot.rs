use crate::cli::PlotArgs;
use crate::trvz::color::pick_colors;
use crate::trvz::waterfall_plot::plot_waterfall;
use crate::trvz::{allele_plot::plot_alleles, input};
use crate::utils::{open_catalog_reader, open_genome_reader, Result};
use pipeplot::generate_image;

pub fn trvz(args: PlotArgs) -> Result<()> {
    let catalog_reader = open_catalog_reader(&args.repeats_path)?;
    let genome_reader = open_genome_reader(&args.genome_path)?;
    let locus = input::get_locus(catalog_reader, genome_reader, &args.tr_id, args.flank_len)?;

    let reads = input::get_reads(&args.reads_path, &locus, args.max_allele_reads)?;
    let colors = pick_colors(&locus.motifs);
    let mut pipe_plot = if args.plot_type == "allele" {
        let allele_seqs = input::get_alleles(&args.bcf_path, &locus)?;
        plot_alleles(&locus, &args.what_to_show, &allele_seqs, &reads, colors)
    } else {
        plot_waterfall(&locus, &args.what_to_show, &reads, &colors)
    };

    if let Some(font_family) = args.font_family {
        pipe_plot.set_font_family(&font_family);
    }

    generate_image(&pipe_plot, &args.output_path)?;
    Ok(())
}
