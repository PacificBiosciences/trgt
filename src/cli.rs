use crate::{
    merge::vcf_writer::OutputType,
    utils::{Genotyper, Result, TrgtPreset, TrgtScoring},
};
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
    let git_describe = env!("VERGEN_GIT_DESCRIBE");
    if git_describe.is_empty() {
        env!("CARGO_PKG_VERSION").to_string()
    } else {
        format!("{}-{}", env!("CARGO_PKG_VERSION"), git_describe)
    }
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
    #[clap(about = "Tandem Repeat VCF Merger")]
    Merge(MergeArgs),
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("merge")))]
#[command(arg_required_else_help(true))]
pub struct MergeArgs {
    #[clap(required = true)]
    #[clap(short = 'v')]
    #[clap(long = "vcf")]
    #[clap(help = "VCF files to merge")]
    #[clap(value_name = "VCF")]
    #[clap(num_args = 1..)]
    #[arg(value_parser = check_file_exists)]
    pub vcfs: Vec<PathBuf>,

    #[clap(short = 'g')]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: Option<PathBuf>,

    #[clap(short = 'o')]
    #[clap(long = "output")]
    #[clap(value_name = "FILE")]
    #[clap(help = "Write output to a file [standard output]")]
    #[arg(value_parser = check_prefix_path)]
    pub output: Option<PathBuf>,

    #[clap(help_heading("Advanced"))]
    #[clap(short = 'O')]
    #[clap(long = "output-type")]
    #[clap(value_name = "OUTPUT_TYPE")]
    #[clap(help = "Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF")]
    #[clap(value_parser = merge_validate_output_type)]
    pub output_type: Option<OutputType>,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "skip-n")]
    #[clap(value_name = "SKIP_N")]
    #[clap(help = "Skip the first N records")]
    pub skip_n: Option<usize>,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "process-n")]
    #[clap(value_name = "process_N")]
    #[clap(help = "Only process N records")]
    pub process_n: Option<usize>,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "print-header")]
    #[clap(help = "Print only the merged header and exit")]
    pub print_header: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "force-single")]
    #[clap(help = "Run even if there is only one file on input")]
    pub force_single: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(hide = true)]
    #[clap(long = "force-samples")]
    #[clap(help = "Resolve duplicate sample names")]
    pub force_samples: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "no-version")]
    #[clap(help = "Do not append version and command line to the header")]
    pub no_version: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(hide = true)]
    #[clap(long = "missing-to-ref")]
    #[clap(help = "Assume genotypes at missing sites are 0/0")]
    pub missing_to_ref: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(hide = true)]
    #[clap(long = "strategy")]
    #[clap(value_name = "STRATEGY")]
    #[clap(help = "Set variant merging strategy to use")]
    #[clap(value_parser(["exact"]))]
    #[clap(default_value = "exact")]
    pub merge_strategy: String,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "quit-on-errors")]
    #[clap(help = "Quit immediately on errors during merging")]
    pub quit_on_error: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "contig")]
    #[clap(value_name = "CONTIG")]
    #[clap(help = "Process only the specified contigs (comma-separated list)")]
    #[clap(value_delimiter = ',')]
    pub contigs: Option<Vec<String>>,
}

#[derive(Parser, Debug, Clone)]
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
    pub output_prefix: PathBuf,

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

    #[clap(long = "preset")]
    #[clap(value_name = "PRESET")]
    #[clap(help = "Parameter preset (wgs or targeted)")]
    #[clap(default_value = "wgs")]
    pub preset: TrgtPreset,

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
    #[clap(default_value_if("preset", "targeted", "cluster"))]
    pub genotyper: Genotyper,

    #[clap(hide = true)]
    #[clap(help_heading("Advanced"))]
    #[clap(long = "aln-scoring")]
    #[clap(value_name = "SCORING")]
    #[clap(
        help = "Scoring function to align to flanks (non-negative values): MATCH,MISM,GAPO,GAPE,KMERLEN,BANDWIDTH"
    )]
    #[clap(default_value = "1,1,5,1,8,6")]
    #[clap(default_value_if("preset", "targeted", "1,1,0,1,1,5000"))]
    #[arg(value_parser = scoring_from_string)]
    pub aln_scoring: TrgtScoring,

    #[clap(hide = true)]
    #[clap(help_heading("Advanced"))]
    #[clap(long = "min-flank-id-frac")]
    #[clap(value_name = "PERC")]
    #[clap(help = "Minimum fraction of matches in a flank sequence to consider it 'found'")]
    #[clap(default_value = "0.7")]
    #[clap(default_value_if("preset", "targeted", "0.8"))]
    #[arg(value_parser = ensure_unit_float)]
    pub min_flank_id_frac: f64,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "flank-len")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(help = "Minimum length of the flanking sequence")]
    #[clap(default_value = "250")]
    #[clap(default_value_if("preset", "targeted", "200"))]
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

    #[clap(hide = true)]
    #[clap(help_heading("Advanced"))]
    #[clap(long = "min-read-quality")]
    #[clap(value_name = "MIN_RQ")]
    #[clap(help = "Minimum HiFi rq value required to use a read for genotyping")]
    #[clap(default_value = "0.98")]
    #[clap(default_value_if("preset", "targeted", "-1.0"))]
    // #[arg(value_parser = ensure_unit_float)]
    pub min_hifi_read_qual: f64,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "disable-bam-output")]
    #[clap(help = "Disable BAM output")]
    pub disable_bam_output: bool,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "max-depth")]
    #[clap(value_name = "MAX_DEPTH")]
    #[clap(help = "Maximum locus depth")]
    #[clap(default_value = "250")]
    #[clap(default_value_if("preset", "targeted", "10000"))]
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
    pub output_path: PathBuf,

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

    #[clap(help_heading("Advanced"))]
    #[clap(long = "max-allele-reads")]
    #[clap(value_name = "MAX_READS")]
    #[clap(help = "Max number of reads per allele to plot")]
    pub max_allele_reads: Option<usize>,
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
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
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

fn check_prefix_path(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if let Some(parent_dir) = path.parent() {
        if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() {
            return Err(format!("Path does not exist: {}", parent_dir.display()));
        }
    }
    Ok(PathBuf::from(s))
}

fn check_image_path(s: &str) -> Result<PathBuf> {
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

fn merge_validate_output_type(s: &str) -> Result<OutputType> {
    let valid_prefixes = ["u", "b", "v", "z"];
    if valid_prefixes.contains(&s) {
        return match s {
            "u" => Ok(OutputType::Bcf {
                is_uncompressed: true,
                level: None,
            }),
            "v" => Ok(OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            }),
            "b" => Ok(OutputType::Bcf {
                is_uncompressed: false,
                level: None,
            }),
            "z" => Ok(OutputType::Vcf {
                is_uncompressed: false,
                level: None,
            }),
            _ => unreachable!(),
        };
    }

    // NOTE: Can't actually set compression level in rust/htslib at the moment
    // if s.len() == 2 {
    //     let (prefix, suffix) = s.split_at(1);
    //     if (prefix == "b" || prefix == "z") && suffix.chars().all(|c| c.is_digit(10)) {
    //         return match prefix {
    //             "b" => Ok(OutputType::Bcf {
    //                 is_uncompressed: false,
    //                 level: Some(suffix.parse().unwrap()),
    //             }),
    //             "z" => Ok(OutputType::Vcf {
    //                 is_uncompressed: false,
    //                 level: Some(suffix.parse().unwrap()),
    //             }),
    //             _ => unreachable!(),
    //         };
    //     } else if (prefix == "u" || prefix == "v") && suffix.chars().all(|c| c.is_digit(10)) {
    //         return Err(format!(
    //             "Error: compression level ({}) cannot be set on uncompressed streams ({})",
    //             suffix, prefix
    //         ));
    //     }
    // }

    Err(format!(
        "Invalid output type: {}. Must be one of u, b, v, z.",
        s
    ))
}
