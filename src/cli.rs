use crate::utils::{Genotyper, Result, TrgtScoring};
use chrono::Datelike;
use clap::{ArgAction, ArgGroup, Parser, Subcommand};
use env_logger::fmt::Color;
use log::{Level, LevelFilter};
use once_cell::sync::Lazy;
use std::{
    io::Write,
    path::{Path, PathBuf},
};

pub static FULL_VERSION: Lazy<String> = Lazy::new(|| {
    format!(
        "{}-{}",
        env!("CARGO_PKG_VERSION"),
        env!("VERGEN_GIT_DESCRIBE")
    )
});

#[derive(Parser)]
#[command(name="trgt",
          author="Egor Dolzhenko <edolzhenko@pacificbiosciences.com>\nGuilherme De Sena Brandine <gbrandine@pacificbiosciences.com>\nTom Mokveld <tmokveld@pacificbiosciences.com>", 
          version=&**FULL_VERSION,
          long_about = None,
          disable_help_subcommand = true,
          after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}\n{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = ArgAction::Count, help = "Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)")]
    pub verbosity: u8,
}

#[derive(Subcommand)]
pub enum Command {
    #[clap(about = "Tandem Repeat Genotyper")]
    Genotype(GenotypeArgs),
    #[clap(about = "Tandem Repeat Plotter")]
    Plot(PlotArgs),
    #[clap(about = "Tandem Repeat Catalog Validator")]
    Validate(ValidateArgs),
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("genotype")))]
#[command(arg_required_else_help(true))]
pub struct GenotypeArgs {
    #[clap(required = true)]
    #[clap(short = 'g')]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reads")]
    #[clap(help = "BAM file with aligned HiFi reads")]
    #[clap(value_name = "READS")]
    #[arg(value_parser = check_file_exists)]
    pub reads_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "repeats")]
    #[clap(help = "BED file with repeat coordinates")]
    #[clap(value_name = "REPEATS")]
    #[arg(value_parser = check_file_exists)]
    pub repeats_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "output-prefix")]
    #[clap(help = "Prefix for output files")]
    #[clap(value_name = "OUTPUT_PREFIX")]
    #[arg(value_parser = check_prefix_path)]
    pub output_prefix: String,

    #[clap(long = "karyotype")]
    #[clap(short = 'k')]
    #[clap(value_name = "KARYOTYPE")]
    #[clap(help = "Sample karyotype (XX or XY or file name)")]
    #[clap(default_value = "XX")]
    pub karyotype: String,

    #[clap(short = 't')]
    #[clap(long = "threads")]
    #[clap(help = "Number of threads")]
    #[clap(value_name = "THREADS")]
    #[clap(default_value = "1")]
    #[arg(value_parser = threads_in_range)]
    pub num_threads: usize,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "sample-name")]
    #[clap(value_name = "SAMPLE_NAME")]
    #[clap(help = "Sample name")]
    #[clap(default_value = None)]
    #[arg(value_parser = check_sample_name_nonempty)]
    pub sample_name: Option<String>,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "genotyper")]
    #[clap(value_name = "GENOTYPER")]
    #[clap(help = "Genotyping algorithm (size or cluster)")]
    #[clap(default_value = "size")]
    pub genotyper: Genotyper,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "aln-scoring")]
    #[clap(value_name = "SCORING")]
    #[clap(
        help = "Scoring function to align to flanks (non-negative values): MATCH,MISM,GAPO,GAPE,KMERLEN,BANDWIDTH"
    )]
    #[clap(default_value = "1,1,5,1,8,6")]
    #[arg(value_parser = scoring_from_string)]
    pub aln_scoring: TrgtScoring,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "min-flank-id-frac")]
    #[clap(value_name = "PERC")]
    #[clap(help = "Minimum fraction of matches in a flank sequence to consider it 'found'")]
    #[clap(default_value = "0.7")]
    #[arg(value_parser = ensure_unit_float)]
    pub min_flank_id_frac: f64,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "flank-len")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(help = "Minimum length of the flanking sequence")]
    #[clap(default_value = "250")]
    pub flank_len: usize,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "output-flank-len")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(help = "Length of flanking sequence to report on output")]
    #[clap(default_value = "50")]
    pub output_flank_len: usize,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "fixed-flanks")]
    #[clap(value_name = "FIXED_FLANKS")]
    #[clap(help = "Keep flank length fixed")]
    pub fixed_flanks: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "min-read-quality")]
    #[clap(value_name = "MIN_RQ")]
    #[clap(help = "Minimum HiFi rq value required to use a read for genotyping")]
    #[clap(default_value = "0.98")]
    // #[arg(value_parser = ensure_unit_float)]
    pub min_hifi_read_qual: f64,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "max-depth")]
    #[clap(value_name = "MAX_DEPTH")]
    #[clap(help = "Maximum locus depth")]
    #[clap(default_value = "250")]
    pub max_depth: usize,
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("plot")))]
#[command(arg_required_else_help(true))]
pub struct PlotArgs {
    #[clap(required = true)]
    #[clap(short = 'g')]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "repeats")]
    #[clap(help = "BED file with repeat coordinates")]
    #[clap(value_name = "REPEATS")]
    #[arg(value_parser = check_file_exists)]
    pub repeats_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'v')]
    #[clap(long = "vcf")]
    #[clap(help = "VCF file generated by TRGT")]
    #[clap(value_name = "VCF")]
    #[arg(value_parser = check_file_exists)]
    pub bcf_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "spanning-reads")]
    #[clap(help = "BAM file with spanning reads generated by TRGT")]
    #[clap(value_name = "SPANNING_READS")]
    #[arg(value_parser = check_file_exists)]
    pub reads_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 't')]
    #[clap(long = "repeat-id")]
    #[clap(help = "ID of the repeat to plot")]
    #[clap(value_name = "REPEAT_ID")]
    pub tr_id: String,

    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "image")]
    #[clap(help = "Output image path")]
    #[clap(value_name = "IMAGE")]
    #[arg(value_parser = check_image_path)]
    pub output_path: String,

    #[clap(help_heading("Plotting"))]
    #[clap(long = "plot-type")]
    #[clap(value_name = "PLOT_TYPE")]
    #[clap(help = "Type of plot to generate")]
    #[clap(value_parser(["allele",  "waterfall"]))]
    #[clap(default_value = "allele")]
    pub plot_type: String,

    #[clap(help_heading("Plotting"))]
    #[clap(long = "show")]
    #[clap(value_name = "SHOW")]
    #[clap(help = "What to show in the plot")]
    #[clap(value_parser(["motifs",  "meth"]))]
    #[clap(default_value = "motifs")]
    pub what_to_show: String,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "flank-len")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(help = "Length of flanking regions")]
    #[clap(default_value = "50")]
    pub flank_len: usize,
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("validate")))]
#[command(arg_required_else_help(true))]
pub struct ValidateArgs {
    #[clap(required = true)]
    #[clap(short = 'g')]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "repeats")]
    #[clap(help = "BED file with repeat coordinates")]
    #[clap(value_name = "REPEATS")]
    #[arg(value_parser = check_file_exists)]
    pub repeats_path: PathBuf,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "flank-len")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(help = "Length of flanking regions")]
    #[clap(default_value = "50")]
    pub flank_len: usize,
}

pub fn init_verbose(args: &Cli) {
    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        _ => LevelFilter::Debug,
    };

    env_logger::Builder::from_default_env()
        .format(|buf, record| {
            let level = record.level();
            let mut style = buf.style();
            match record.level() {
                Level::Error => style.set_color(Color::Red),
                Level::Warn => style.set_color(Color::Yellow),
                Level::Info => style.set_color(Color::Green),
                Level::Debug => style.set_color(Color::Blue),
                Level::Trace => style.set_color(Color::Cyan),
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                style.value(level),
                record.args()
            )
        })
        .filter_level(filter_level)
        .init();
}

fn check_prefix_path(s: &str) -> Result<String> {
    let path = Path::new(s);
    if let Some(parent_dir) = path.parent() {
        if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() {
            return Err(format!("Path does not exist: {}", parent_dir.display()));
        }
    }
    Ok(s.to_string())
}

fn check_image_path(s: &str) -> Result<String> {
    let prefix_check = check_prefix_path(s)?;
    let path = Path::new(s);
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("svg") | Some("png") | Some("pdf") => Ok(prefix_check),
        _ => Err("Image must have an extension of .svg, .png, or .pdf".to_string()),
    }
}

fn threads_in_range(s: &str) -> Result<usize> {
    let thread: usize = s
        .parse()
        .map_err(|_| format!("`{}` is not a valid thread number", s))?;
    if thread >= 1 {
        Ok(thread)
    } else {
        Err("Number of threads must be at least 1".into())
    }
}

fn check_file_exists(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if !path.exists() {
        Err(format!("File does not exist: {}", path.display()))
    } else {
        Ok(path.to_path_buf())
    }
}

fn check_sample_name_nonempty(s: &str) -> Result<String> {
    if s.trim().is_empty() {
        Err("Sample name cannot be an empty string".to_string())
    } else {
        Ok(s.to_string())
    }
}

fn ensure_unit_float(s: &str) -> Result<f64> {
    let value = s
        .parse::<f64>()
        .map_err(|e| format!("Could not parse float: {}", e))?;
    if !(0.0..=1.0).contains(&value) {
        Err(format!(
            "The value must be between 0.0 and 1.0, got: {}",
            value
        ))
    } else {
        Ok(value)
    }
}

fn scoring_from_string(s: &str) -> Result<TrgtScoring> {
    const NUM_EXPECTED_VALUES: usize = 6;
    let values: Vec<i32> = s.split(',').filter_map(|x| x.parse().ok()).collect();
    if values.len() != NUM_EXPECTED_VALUES {
        return Err(format!(
            "Expected {} comma-separated values in scoring. Got {} -> {}",
            NUM_EXPECTED_VALUES,
            values.len(),
            s
        ));
    }

    if values.iter().any(|&val| val < 0) {
        return Err(format!(
            "Negative values are not allowed in scoring. Got {}.",
            s
        ));
    }

    Ok(TrgtScoring {
        match_scr: values[0],
        mism_scr: values[1],
        gapo_scr: values[2],
        gape_scr: values[3],
        kmer_len: values[4] as usize,
        bandwidth: values[5] as usize,
    })
}
