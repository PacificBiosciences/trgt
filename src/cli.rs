use crate::{
    merge::vcf_reader::VcfReadMode,
    preflight_fields,
    utils::{
        BamSources, Genotyper, InputSource, Karyotype, OutputType, Result, TrgtPreset, TrgtScoring,
    },
};
use chrono::Datelike;
use clap::{builder::ArgPredicate, ArgAction, ArgGroup, Parser, Subcommand, ValueEnum};
use log::{Level, LevelFilter};
use owo_colors::{
    colors::{Blue, Green, Magenta, Red, Yellow},
    OwoColorize, Stream, Style,
};
use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

#[cfg(has_git_describe)]
pub const FULL_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "-", env!("VERGEN_GIT_DESCRIBE"));

#[cfg(not(has_git_describe))]
pub const FULL_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[command(name="trgt",
          author="Egor Dolzhenko <edolzhenko@pacificbiosciences.com>\nGuilherme De Sena Brandine <gbrandine@pacificbiosciences.com>\nTom Mokveld <tmokveld@pacificbiosciences.com>", 
          version=FULL_VERSION,
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
    /// Enable or disable color output in logging
    #[arg(long, value_enum, default_value_t = Color::Auto, global = true, help_heading = "Advanced")]
    color: Color,

    /// Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
    #[arg(
        short = 'v',
        long = "verbose",
        action = ArgAction::Count,
        global = true
    )]
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
    #[clap(about = "Locus-specific Deep Dive")]
    Deepdive(DeepdiveArgs),
}

impl Command {
    pub fn name(&self) -> &'static str {
        match self {
            Command::Genotype(_) => "genotype",
            Command::Plot(_) => "plot",
            Command::Validate(_) => "validate",
            Command::Merge(_) => "merge",
            Command::Deepdive(_) => "deepdive",
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(group(ArgGroup::new("genotype")))]
#[command(arg_required_else_help(true))]
pub struct GenotypeArgs {
    /// Path to reference genome FASTA
    #[arg(short = 'g', long = "genome", value_name = "FASTA", required = true)]
    pub genome_src: InputSource,

    /// BAM file with aligned HiFi reads
    #[arg(
        short = 'r',
        long = "reads",
        value_name = "READS",
        required = true,
        num_args = 1
    )]
    pub reads_src: InputSource,

    /// Optional BAM file with fail reads, also disables the default rq filter unless explicitly set with --min-read-quality
    #[arg(long = "fail-reads", value_name = "FAIL_READS", num_args = 1)]
    pub fail_reads_src: Option<InputSource>,

    /// BED file with repeat coordinates
    #[arg(short = 'b', long = "repeats", value_name = "REPEATS", required = true)]
    pub repeats_src: InputSource,

    /// Prefix for output files (.vcf/.bcf/.vcf.gz/.bcf.gz and .spanning.bam)
    #[arg(
        short = 'o',
        long = "output-prefix",
        value_name = "OUTPUT_PREFIX",
        value_parser = check_prefix_path,
        required = true
    )]
    pub output_prefix: PathBuf,

    /// Sample karyotype (XX or XY or file name)
    #[arg(
        short = 'k',
        long = "karyotype",
        value_name = "KARYOTYPE",
        default_value = "XX",
        value_parser = Karyotype::new
    )]
    pub karyotype: Karyotype,

    /// Number of threads
    #[arg(
        short = 't',
        long = "threads",
        value_name = "THREADS",
        default_value = "1",
        value_parser = threads_in_range
    )]
    pub num_threads: usize,

    /// Parameter preset (wgs or targeted)
    #[arg(long = "preset", value_name = "PRESET", default_value = "wgs")]
    pub preset: TrgtPreset,

    /// Sample name
    #[arg(
        long = "sample-name",
        value_name = "SAMPLE_NAME",
        default_value = None,
        value_parser = check_sample_name_nonempty,
        help_heading = "Advanced"
    )]
    pub sample_name: Option<String>,

    /// Genotyping algorithm (size or cluster)
    #[arg(
        long = "genotyper",
        value_name = "GENOTYPER",
        default_value = "size",
        default_value_if("preset", "targeted", Some("cluster")),
        help_heading = "Advanced"
    )]
    pub genotyper: Genotyper,

    /// Scoring function to align to flanks (non-negative values): MISM,GAPO,GAPE
    #[arg(
         long = "aln-scoring",
         value_name = "SCORING",
         default_value = "2,5,1",
         default_value_if("preset", "targeted", Some("1,0,1")),
         value_parser = scoring_from_string,
         help_heading = "Advanced",
         hide = true
     )]
    pub aln_scoring: TrgtScoring,

    /// Minimum fraction of matches in a flank sequence to consider it 'found'
    #[arg(
        long = "min-flank-id-frac",
        value_name = "PERC",
        default_value = "0.7",
        default_value_if("preset", "targeted", Some("0.8")),
        value_parser = ensure_unit_float,
        help_heading = "Advanced",
        hide = true
    )]
    pub min_flank_id_frac: f64,

    /// Minimum length of the flanking sequence
    #[arg(
        long = "flank-len",
        value_name = "FLANK_LEN",
        default_value = "250",
        default_value_if("preset", "targeted", "200"),
        help_heading = "Advanced"
    )]
    pub flank_len: usize,

    /// Length of flanking sequence to report on output
    #[arg(
        long = "output-flank-len",
        value_name = "FLANK_LEN",
        default_value = "50",
        help_heading = "Advanced IO"
    )]
    pub output_flank_len: usize,

    /// Keep flank length fixed
    #[arg(
        long = "fixed-flanks",
        value_name = "FIXED_FLANKS",
        help_heading = "Advanced",
        hide = true
    )]
    pub fixed_flanks: bool,

    /// Minimum HiFi rq value required to use a read for genotyping
    #[arg(
        long = "min-read-quality",
        value_name = "MIN_RQ",
        default_value = "0.98",
        default_value_if("preset", "targeted", Some("-1.0")),
        default_value_if("fail_reads_src", ArgPredicate::IsPresent, Some("-1.0")),
        help_heading = "Advanced",
        hide = true
    )]
    pub min_hifi_read_qual: f32,

    /// Disable BAM output
    #[arg(long = "disable-bam-output", help_heading = "Advanced IO")]
    pub disable_bam_output: bool,

    /// Maximum locus depth
    #[arg(
        long = "max-depth",
        value_name = "MAX_DEPTH",
        default_value = "250",
        default_value_if("preset", "targeted", Some("10000")),
        help_heading = "Advanced"
    )]
    pub max_depth: usize,

    /// Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "OUTPUT_TYPE",
        value_parser = validate_output_type,
        help_heading = "Advanced IO"
    )]
    pub output_type: Option<OutputType>,

    /// Maximum span of a locus group in base pairs. Loci separated by more than this distance will be in different groups
    #[arg(
        long = "max-group-span",
        value_name = "MAX_GROUP_SPAN",
        default_value = "50000",
        help_heading = "Advanced IO",
        hide = true
    )]
    pub max_group_span: u32,

    /// Maximum size of a locus group
    #[arg(
        long = "max-group-size",
        value_name = "MAX_GROUP_SIZE",
        default_value = "10000",
        help_heading = "Advanced IO",
        hide = true
    )]
    pub max_group_size: usize,

    /// Number of threads for querying input BAM files. Defaults to half the number of analysis threads, with a max of 8
    #[arg(
        long = "fetcher-threads",
        value_name = "THREADS",
        value_parser = threads_in_range,
        help_heading = "Advanced IO"
    )]
    pub fetcher_threads: Option<usize>,

    /// Number of threads for decompressing input files
    #[arg(
        long = "decompression-threads",
        value_name = "THREADS",
        default_value = "1",
        value_parser = threads_in_range,
        help_heading = "Advanced IO",
        hide = true
    )]
    pub decompression_threads: usize,

    /// Number of threads for BAM compression
    #[arg(
        long = "bam-compression-threads",
        value_name = "THREADS",
        value_parser = threads_in_range,
        help_heading = "Advanced IO",
        hide = true
    )]
    pub bam_compression_threads: Option<usize>,

    /// Buffer size for the locus group channel
    #[arg(
        long = "group-channel-buffer",
        value_name = "SIZE",
        help_heading = "Advanced IO",
        hide = true
    )]
    pub group_channel_buffer_size: Option<usize>,

    /// Buffer size for the locus-with-reads channel
    #[arg(
        long = "locus-channel-buffer",
        value_name = "SIZE",
        default_value = "4096",
        help_heading = "Advanced IO",
        hide = true
    )]
    pub locus_channel_buffer_size: usize,

    /// Buffer size for the results channel
    #[arg(
        long = "result-channel-buffer",
        value_name = "SIZE",
        default_value = "2048",
        help_heading = "Advanced IO",
        hide = true
    )]
    pub result_channel_buffer_size: usize,

    /// Skip phasing annotations
    #[arg(long = "skip-phase-annotation", help_heading = "Advanced", hide = true)]
    pub skip_phase_annotation: bool,
}

impl GenotypeArgs {
    pub fn read_sources(&self) -> BamSources {
        BamSources::new(self.reads_src.clone(), self.fail_reads_src.clone())
    }

    pub fn preflight(&self) -> Result<()> {
        preflight_fields!(self, genome_src, reads_src, fail_reads_src, repeats_src)
    }
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("plot")))]
#[command(arg_required_else_help(true))]
pub struct PlotArgs {
    /// Path to reference genome FASTA
    #[arg(short = 'g', long = "genome", value_name = "FASTA", required = true)]
    pub genome_src: InputSource,

    /// BED file with repeat coordinates
    #[arg(short = 'b', long = "repeats", value_name = "REPEATS", required = true)]
    pub repeats_src: InputSource,

    /// VCF file generated by TRGT
    #[arg(short = 'f', long = "vcf", value_name = "VCF", required = true)]
    pub bcf_src: InputSource,

    /// BAM file with spanning reads generated by TRGT
    #[arg(
        short = 'r',
        long = "spanning-reads",
        value_name = "SPANNING_READS",
        required = true
    )]
    pub reads_src: InputSource,

    /// ID of the repeat to plot
    #[arg(
        short = 'i',
        long = "repeat-id",
        value_name = "REPEAT_ID",
        required = true
    )]
    pub tr_id: String,

    /// Output image path
    #[arg(
        short = 'o',
        long = "image",
        value_name = "IMAGE",
        value_parser = check_image_path,
        required = true
    )]
    pub output_path: PathBuf,

    /// Type of plot to generate
    #[arg(
        long = "plot-type",
        value_name = "PLOT_TYPE",
        value_parser = ["allele", "waterfall"],
        default_value = "allele",
        help_heading = "Plotting"
    )]
    pub plot_type: String,

    /// What to show in the plot
    #[arg(
        long = "show",
        value_name = "SHOW",
        value_parser = ["motifs", "meth"],
        default_value = "motifs",
        help_heading = "Plotting"
    )]
    pub what_to_show: String,

    /// Horizontally compress the plot
    #[arg(long = "squished", help_heading = "Plotting")]
    pub is_squished: bool,

    /// Font family to use for text elements (default: Roboto Mono)
    #[arg(long = "font-family", value_name = "FONT", help_heading = "Plotting")]
    pub font_family: Option<String>,

    /// Length of flanking regions
    #[arg(
        long = "flank-len",
        value_name = "FLANK_LEN",
        default_value = "50",
        help_heading = "Advanced"
    )]
    pub flank_len: usize,

    /// Max number of reads per allele to plot
    #[arg(
        long = "max-allele-reads",
        value_name = "MAX_READS",
        help_heading = "Advanced"
    )]
    pub max_allele_reads: Option<usize>,
}

impl PlotArgs {
    pub fn preflight(&self) -> Result<()> {
        preflight_fields!(self, genome_src, repeats_src, bcf_src, reads_src)
    }
}

#[derive(Parser, Debug)]
#[command(group(
    ArgGroup::new("input")
        .required(true)
        .args(["vcfs", "vcf_list"]),
))]
#[command(arg_required_else_help(true))]
pub struct MergeArgs {
    /// VCF files to merge
    #[arg(
        long = "vcf",
        value_name = "VCF",
        num_args = 1..,
        value_parser = check_file_exists
    )]
    pub vcfs: Option<Vec<PathBuf>>,

    /// File containing paths of VCF files to merge (one per line)
    #[arg(
        long = "vcf-list",
        value_name = "VCF_LIST",
        value_parser = check_file_exists
    )]
    pub vcf_list: Option<PathBuf>,

    /// Path to reference genome FASTA
    #[arg(short = 'g', long = "genome", value_name = "FASTA")]
    pub genome_src: Option<InputSource>,

    /// Write output to a file [standard output]
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        value_parser = check_prefix_path
    )]
    pub output: Option<PathBuf>,

    /// Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "OUTPUT_TYPE",
        value_parser = validate_output_type,
        help_heading = "Advanced"
    )]
    pub output_type: Option<OutputType>,

    /// Skip the first N records
    #[arg(long = "skip-n", value_name = "SKIP_N", help_heading = "Advanced")]
    pub skip_n: Option<usize>,

    /// Only process N records
    #[arg(
        long = "process-n",
        value_name = "PROCESS_N",
        help_heading = "Advanced"
    )]
    pub process_n: Option<usize>,

    /// Print only the merged header and exit
    #[arg(long = "print-header", help_heading = "Advanced")]
    pub print_header: bool,

    /// Run even if there is only one file on input
    #[arg(long = "force-single", help_heading = "Advanced")]
    pub force_single: bool,

    /// Do not append version and command line to the header
    #[arg(long = "no-version", help_heading = "Advanced")]
    pub no_version: bool,

    /// Quit immediately on errors during merging
    #[arg(long = "quit-on-errors", help_heading = "Advanced")]
    pub quit_on_error: bool,

    /// Process only the specified contigs (comma-separated list)
    #[arg(
        long = "contig",
        value_name = "CONTIG",
        value_delimiter = ',',
        help_heading = "Advanced"
    )]
    pub contigs: Option<Vec<String>>,

    /// Number of threads for (de)compressing input/output VCF files (a threadpool is shared between all readers and the writer)
    #[arg(
        short = 't',
        long = "threads",
        value_name = "THREADS",
        default_value = "2",
        value_parser = threads_in_range,
        help_heading = "Advanced"
    )]
    pub threads: usize,

    /// Stream VCFs without loading their indexes (contig order must match across inputs)
    #[arg(
        long = "no-index",
        help_heading = "Advanced",
        value_enum,
        default_value_t = VcfReadMode::Indexed,
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "streaming",
        value_name = "MODE"
    )]
    pub no_index: VcfReadMode,

    /// Write index for the output compressed VCF/BCF file
    #[arg(short = 'W', long = "write-index", help_heading = "Advanced")]
    pub write_index: bool,
}

impl MergeArgs {
    pub fn process_vcf_paths(&self) -> Result<Vec<PathBuf>> {
        match (&self.vcfs, &self.vcf_list) {
            (Some(vcfs), None) => Ok(vcfs.clone()),
            (None, Some(list_path)) => Self::read_vcf_paths_from_file(list_path),
            _ => unreachable!("Either --vcf or --vcf-list is provided, never both"),
        }
    }

    fn read_vcf_paths_from_file(path: &Path) -> Result<Vec<PathBuf>> {
        let file = File::open(path)
            .map_err(|e| format!("Failed to open VCF list file {}: {}", path.display(), e))?;
        let reader = BufReader::new(file);

        let mut paths = Vec::new();
        for (line_num, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| format!("Error reading line {}: {}", line_num + 1, e))?;
            let trimmed = line.trim();
            // Skip empty or comment lines
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let path = PathBuf::from(trimmed);
            if !path.exists() {
                Err(format!("VCF file does not exist: {}", path.display()))?;
            }
            paths.push(path);
        }

        if paths.is_empty() {
            Err("No VCF paths found in the input file".to_string())?;
        }

        Ok(paths)
    }

    pub fn preflight(&self) -> Result<()> {
        preflight_fields!(self, genome_src)
    }
}

#[derive(Parser, Debug, Clone)]
#[command(group(ArgGroup::new("deepdive")))]
#[command(arg_required_else_help(true))]
pub struct DeepdiveArgs {
    /// Path to reference genome FASTA
    #[arg(short = 'g', long = "genome", value_name = "FASTA", required = true)]
    pub genome_src: InputSource,

    /// VCF file generated by trgt genotype
    #[arg(short = 'f', long = "vcf", value_name = "VCF", required = true)]
    pub bcf_src: InputSource,

    /// BED file with repeat coordinates
    #[arg(short = 'b', long = "repeats", value_name = "REPEATS", required = true)]
    pub repeats_src: InputSource,

    /// BAM file with spanning reads generated by trgt genotype
    #[arg(
        short = 'r',
        long = "spanning-reads",
        value_name = "SPANNING_READS",
        required = true
    )]
    pub reads_src: InputSource,

    /// ID of the repeat to realign
    #[arg(
        short = 'i',
        long = "repeat-id",
        value_name = "REPEAT_ID",
        required = true
    )]
    pub tr_id: String,

    /// Prefix for output files (.fasta, .bed, and .bam)
    #[arg(
        short = 'o',
        long = "output-prefix",
        value_name = "OUTPUT_PREFIX",
        value_parser = check_prefix_path,
        required = true
    )]
    pub output_prefix: PathBuf,
}

impl DeepdiveArgs {
    pub fn preflight(&self) -> Result<()> {
        preflight_fields!(self, genome_src, bcf_src, repeats_src, reads_src)
    }
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("validate")))]
#[command(arg_required_else_help(true))]
pub struct ValidateArgs {
    /// Path to reference genome FASTA
    #[arg(short = 'g', long = "genome", value_name = "FASTA", required = true)]
    pub genome_src: InputSource,

    /// BED file with repeat coordinates
    #[arg(short = 'b', long = "repeats", value_name = "REPEATS", required = true)]
    pub repeats_src: InputSource,

    /// Length of flanking regions
    #[arg(
        long = "flank-len",
        value_name = "FLANK_LEN",
        default_value = "50",
        help_heading = "Advanced"
    )]
    pub flank_len: usize,
}

impl ValidateArgs {
    pub fn preflight(&self) -> Result<()> {
        preflight_fields!(self, genome_src, repeats_src)
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Color {
    Always,
    Auto,
    Never,
}

impl Color {
    fn apply(self) {
        match self {
            Color::Always => owo_colors::set_override(true),
            Color::Auto => {}
            Color::Never => owo_colors::set_override(false),
        }
    }
}

pub fn init_verbose(args: &Cli) {
    args.color.apply();

    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::Builder::from_default_env()
        .format(format_log)
        .filter_level(filter_level)
        .init();
}

#[inline(always)]
fn level_style(level: Level) -> (&'static str, Style) {
    match level {
        Level::Error => ("ERROR", Style::new().fg::<Red>().bold()),
        Level::Warn => ("WARN", Style::new().fg::<Yellow>()),
        Level::Info => ("INFO", Style::new().fg::<Green>()),
        Level::Debug => ("DEBUG", Style::new().fg::<Blue>()),
        Level::Trace => ("TRACE", Style::new().fg::<Magenta>()),
    }
}

fn format_log(buf: &mut env_logger::fmt::Formatter, record: &log::Record) -> std::io::Result<()> {
    let (label, style) = level_style(record.level());
    let ts = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
    let painted_label = label.if_supports_color(Stream::Stderr, |t| style.style(t));
    writeln!(buf, "{ts} [{}] - {}", painted_label, record.args())
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
    const NUM_EXPECTED_VALUES: usize = 3;
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
        mism_scr: values[0],
        gapo_scr: values[1],
        gape_scr: values[2],
    })
}

fn validate_output_type(s: &str) -> Result<OutputType> {
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
