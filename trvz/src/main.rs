use align::{align_reads, align_to_flanks};
use cli::{get_args, handle_error_and_exit};
use flate2::read::GzDecoder;
use genotype_plot::plot_genotype;
use pipe_plot::{DisplayParams, PipePlot};
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::Read as ioRead;
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;
use waterfall_plot::plot_waterfall;
mod align;
mod cli;
mod genotype_plot;
mod hmm;
mod hmm_defs;
mod input;
mod label_hmm;
mod label_motifs;
mod locus;
mod pipe_plot;
mod png;
mod read;
mod struc;
mod svg;
mod utils;
mod waterfall_plot;

pub type Result<T> = std::result::Result<T, String>;

enum FileType {
    Svg,
    Png,
    Pdf,
    Unknown,
}

impl FileType {
    fn from_path(path: &str) -> Self {
        if path.ends_with(".svg") {
            FileType::Svg
        } else if path.ends_with(".png") {
            FileType::Png
        } else if path.ends_with(".pdf") {
            FileType::Pdf
        } else {
            FileType::Unknown
        }
    }
}

fn create_image(plot: PipePlot, path: String) -> io::Result<()> {
    let file_type = FileType::from_path(&path);

    match file_type {
        FileType::Unknown => {
            handle_error_and_exit(format!(
                "Output path {} must end with .svg, .png, or .pdf",
                path
            ));
        }
        _ => (),
    }

    let temp_dir = PathBuf::from(&path);
    let temp_dir = temp_dir.parent().unwrap();
    let svg_temp_path = NamedTempFile::new_in(temp_dir).unwrap().into_temp_path();
    svg::generate(&plot, svg_temp_path.to_str().unwrap());

    match file_type {
        FileType::Svg => {
            svg_temp_path.persist(path).unwrap();
        }
        FileType::Png => {
            png::render(svg_temp_path.to_str().unwrap(), &path);
        }
        FileType::Pdf => {
            let svg = std::fs::read_to_string(svg_temp_path.to_str().unwrap()).unwrap();
            let mut options = usvg::Options::default();
            options.fontdb = usvg::fontdb::Database::new();
            options.fontdb.load_system_fonts();

            let svg = usvg::Tree::from_str(&svg, &options.to_ref()).expect("Error parsing SVG");
            let pdf = svg2pdf::convert_tree(&svg, svg2pdf::Options::default());
            std::fs::write(path, pdf).unwrap();
        }
        FileType::Unknown => unreachable!(),
    }
    Ok(())
}

fn open_catalog_reader(path: &PathBuf) -> Result<BufReader<Box<dyn ioRead>>> {
    fn get_format(path: &Path) -> Option<&'static str> {
        let path_str = path.to_string_lossy();
        let formats = ["bed", "bed.gz", "bed.gzip"];
        formats
            .iter()
            .find(|&&format| path_str.ends_with(format))
            .copied()
    }
    let file = File::open(path).map_err(|e| e.to_string())?;
    match get_format(path) {
        Some("bed.gz") | Some("bed.gzip") => {
            let gz_decoder = GzDecoder::new(file);
            if gz_decoder.header().is_some() {
                Ok(BufReader::new(Box::new(gz_decoder)))
            } else {
                Err(format!("Invalid gzip header: {}", path.to_string_lossy()).into())
            }
        }
        Some("bed") => Ok(BufReader::new(Box::new(file))),
        _ => Err(format!(
            "Unknown bed format: {}. Supported formats are: .bed or .bed.gz(ip)",
            path.to_string_lossy()
        )
        .into()),
    }
}

fn main() -> Result<()> {
    let cli_params = get_args();
    log::info!(
        "Running {}-{}",
        env!("CARGO_PKG_NAME"),
        *crate::cli::FULL_VERSION
    );

    let display_params = DisplayParams {
        what_to_show: cli_params.what_to_show,
        max_width: 750,
    };

    let catalog_reader = open_catalog_reader(&cli_params.repeats_path)
        .unwrap_or_else(|err| handle_error_and_exit(err));

    let locus = input::get_locus(
        &cli_params.genome_path,
        catalog_reader,
        &cli_params.tr_id,
        cli_params.flank_len,
    )
    .unwrap_or_else(|err| handle_error_and_exit(err.into()));

    let reads = input::get_reads(&cli_params.reads_path, &locus)
        .unwrap_or_else(|err| handle_error_and_exit(err.into()));
    let pipe_plot = if cli_params.plot_type == "allele" {
        let genotype = input::get_genotype(&cli_params.bcf_path, &locus)
            .unwrap_or_else(|err| handle_error_and_exit(err.into()));
        let aligns = align_reads(&genotype, reads);
        plot_genotype(&locus, &display_params, &genotype, &aligns)
    } else {
        let aligns = align_to_flanks(&locus, reads);
        plot_waterfall(&locus, &display_params, &aligns)
    };

    create_image(pipe_plot, cli_params.output_path)
        .unwrap_or_else(|err| handle_error_and_exit(err.to_string()));
    log::info!("{} end", env!("CARGO_PKG_NAME"));
    Ok(())
}
