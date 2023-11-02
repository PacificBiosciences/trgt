use crate::locate::TrgtScoring;
use crate::locus::Genotyper;
use chrono::Datelike;
use clap::Parser;
use env_logger::fmt::Color;
use lazy_static::lazy_static;
use log::{Level, LevelFilter};
use std::io::Write;
use std::path::{Path, PathBuf};

lazy_static! {
    pub static ref FULL_VERSION: String = format!(
        "{}-{}",
        env!("CARGO_PKG_VERSION"),
        env!("VERGEN_GIT_DESCRIBE")
    );
}

#[derive(Parser)]
#[command(name="trgt",
          author="Egor Dolzhenko <edolzhenko@pacificbiosciences.com>\nGuilherme De Sena Brandine <gbrandine@pacificbiosciences.com>\nTom Mokveld <tmokveld@pacificbiosciences.com>", 
          version=&**FULL_VERSION,
          about="Tandem Repeat Genotyper", 
          long_about = None,
          after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
          This program comes with ABSOLUTELY NO WARRANTY; it is intended for
          Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
#[command(arg_required_else_help(true))]
pub struct CliParams {
    #[clap(required = true)]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: PathBuf,

    #[clap(required = true)]
    #[clap(long = "reads")]
    #[clap(help = "BAM file with aligned HiFi reads")]
    #[clap(value_name = "READS")]
    #[arg(value_parser = check_file_exists)]
    pub reads_path: PathBuf,

    #[clap(required = true)]
    #[clap(long = "repeats")]
    #[clap(help = "BED file with repeat coordinates")]
    #[clap(value_name = "REPEATS")]
    #[arg(value_parser = check_file_exists)]
    pub repeats_path: PathBuf,

    #[clap(required = true)]
    #[clap(long = "output-prefix")]
    #[clap(help = "Prefix for output files")]
    #[clap(value_name = "OUTPUT_PREFIX")]
    #[arg(value_parser = check_prefix_path)]
    pub output_prefix: String,

    #[clap(long = "threads")]
    #[clap(help = "Number of threads")]
    #[clap(value_name = "THREADS")]
    #[clap(default_value = "1")]
    #[arg(value_parser = threads_in_range)]
    pub num_threads: usize,

    #[clap(long = "karyotype")]
    #[clap(value_name = "KARYOTYPE")]
    #[clap(help = "Sample karyotype (XX or XY or file name)")]
    #[clap(default_value = "XX")]
    pub karyotype: String,

    #[clap(long = "max-depth")]
    #[clap(value_name = "MAX_DEPTH")]
    #[clap(help = "Maximum locus depth")]
    #[clap(default_value = "250")]
    pub max_depth: usize,

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
    #[clap(long = "min-flank-id-perc")]
    #[clap(value_name = "PERC")]
    #[clap(help = "Minimum fraction of matches in a flank sequence to consider it 'found'")]
    #[clap(default_value = "0.7")]
    pub min_flank_id_frac: f32,

    #[clap(help_heading("Advanced"))]
    #[clap(long = "flank-len")]
    #[clap(value_name = "FLANK_LEN")]
    #[clap(help = "Minimum length of the flanking sequence")]
    #[clap(default_value = "250")]
    pub flank_len: usize,

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
    pub min_hifi_read_qual: f64,

    #[clap(help_heading("Advanced"))]
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,
}

pub fn get_cli_params() -> CliParams {
    let args = CliParams::parse();
    check_args(args)
}

fn check_args(args: CliParams) -> CliParams {
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

    args
}

pub fn handle_error_and_exit(err: String) -> ! {
    log::error!("{}", err);
    std::process::exit(1);
}

fn check_prefix_path(s: &str) -> Result<String, String> {
    let path = Path::new(s);
    if let Some(parent_dir) = path.parent() {
        if !parent_dir.as_os_str().is_empty() && !parent_dir.exists() {
            return Err(format!("Path does not exist: {}", parent_dir.display()));
        }
    }
    Ok(s.to_string())
}

fn threads_in_range(s: &str) -> Result<usize, String> {
    let thread: usize = s
        .parse()
        .map_err(|_| format!("`{}` is not a valid thread number", s))?;
    if thread >= 1 {
        Ok(thread)
    } else {
        Err("Number of threads must be at least 1".into())
    }
}

fn check_file_exists(s: &str) -> Result<PathBuf, String> {
    let path = Path::new(s);
    if !path.exists() {
        Err(format!("File does not exist: {}", path.display()))
    } else {
        Ok(path.to_path_buf())
    }
}

fn scoring_from_string(s: &str) -> Result<TrgtScoring, String> {
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
