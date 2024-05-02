use crate::cli::PlotArgs;
use crate::trvz::{
    align::{align_reads, align_to_flanks},
    genotype_plot::plot_genotype,
    input, pdf,
    pipe_plot::{DisplayParams, PipePlot},
    png, svg,
    waterfall_plot::plot_waterfall,
};
use crate::utils::{open_catalog_reader, open_genome_reader, Result};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

#[derive(Debug, PartialEq)]
enum FileType {
    Svg,
    Png,
    Pdf,
}

impl FileType {
    fn from_extension(extension: &str) -> Option<Self> {
        match extension.to_lowercase().as_str() {
            "svg" => Some(FileType::Svg),
            "png" => Some(FileType::Png),
            "pdf" => Some(FileType::Pdf),
            _ => None,
        }
    }
}

fn create_image(plot: PipePlot, path: String) -> Result<()> {
    let extension = Path::new(&path)
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or("Failed to get file extension")?;
    let file_type = FileType::from_extension(extension).ok_or("Unsupported file extension")?;

    let temp_dir = PathBuf::from(&path);
    let temp_dir = temp_dir.parent().unwrap();
    let svg_temp_path = NamedTempFile::new_in(temp_dir).unwrap().into_temp_path();

    svg::generate(&plot, svg_temp_path.to_str().unwrap());
    match file_type {
        FileType::Svg => svg_temp_path.persist(path).unwrap(),
        FileType::Png => png::render(svg_temp_path.to_str().unwrap(), &path),
        FileType::Pdf => pdf::render(svg_temp_path.to_str().unwrap(), &path),
    }
    Ok(())
}

pub fn trvz(args: PlotArgs) -> Result<()> {
    let display_params = DisplayParams {
        what_to_show: args.what_to_show,
        max_width: 750,
    };

    let catalog_reader = open_catalog_reader(&args.repeats_path)?;
    let genome_reader = open_genome_reader(&args.genome_path)?;

    let locus = input::get_locus(catalog_reader, genome_reader, &args.tr_id, args.flank_len)?;

    let reads = input::get_reads(&args.reads_path, &locus)?;
    let pipe_plot = if args.plot_type == "allele" {
        let genotype = input::get_genotype(&args.bcf_path, &locus)?;
        let aligns = align_reads(&genotype, reads);
        plot_genotype(&locus, &display_params, &genotype, &aligns)
    } else {
        let aligns = align_to_flanks(&locus, reads);
        plot_waterfall(&locus, &display_params, &aligns)
    };

    create_image(pipe_plot, args.output_path)?;
    Ok(())
}
