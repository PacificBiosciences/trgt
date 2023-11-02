//! # Tandem Repeat Genotyping Tool (TRGT)
//! TRGT is a program for analysis of tandem repeats designed for [PacBio HiFi
//! sequencing data](https://www.pacb.com/technology/hifi-sequencing/).
//!
//! In order to use TRGT, you will need (a) a reference genome in the FASTA
//! format, (b) a set of repeat definitions in the BED format, and (c) a BAM
//! file containing aligned reads. TRGT will output two files: (a) a VCF file
//! containing repeat allele sequences and their mean methylation levels and (b)
//! a BAM file with pieces of reads overlapping each repeat allele.
//!
//! TRGT can be run like so:
//! ```bash
//!  ./trgt --genome reference.fasta \
//!         --repeats repeats.bed \
//!         --reads read_alignments.bam \
//!         --output-prefix sample
//! ```

use crate::read_output::BamWriter;
use crate::vcf::VcfWriter;
use cli::{get_cli_params, handle_error_and_exit};
use flate2::read::GzDecoder;
use locus::Karyotype;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufReader, Read as ioRead};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::{thread, time};
use threadpool::ThreadPool;
use workflows::analyze_tr;
mod cli;
mod cluster;
mod genotype;
mod label;
mod locate;
mod locus;
mod read_output;
mod reads;
mod snp;
mod utils;
mod vcf;
mod workflows;

pub type Result<T> = std::result::Result<T, String>;

struct ThreadLocalData {
    bam: RefCell<Option<bam::IndexedReader>>,
}

thread_local! {
    static LOCAL: ThreadLocalData = ThreadLocalData {
        bam: RefCell::new(None),
    };
}

pub fn get_bam_header(bam_path: &PathBuf) -> Result<bam::Header> {
    let bam = match bam::IndexedReader::from_path(bam_path) {
        Ok(reader) => reader,
        Err(e) => return Err(format!("Failed to create bam reader: {}", e)),
    };
    let bam_header = bam::Header::from_template(bam.header());
    Ok(bam_header)
}

fn get_sample_name(reads_path: &Path) -> String {
    reads_path
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string()
}

fn create_writer<T, F>(output_prefix: &str, output_suffix: &str, f: F) -> Result<T>
where
    F: FnOnce(&str) -> Result<T>,
{
    let output_path = format!("{}.{}", output_prefix, output_suffix);
    f(&output_path).map_err(|e| {
        eprintln!("Error creating writer: {}", e);
        e
    })
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
                Err(format!("Invalid gzip header: {}", path.to_string_lossy()))
            }
        }
        Some("bed") => Ok(BufReader::new(Box::new(file))),
        _ => Err(format!(
            "Unknown bed format: {}. Supported formats are: .bed or .bed.gz(ip)",
            path.to_string_lossy()
        )),
    }
}

fn open_genome_reader(path: &PathBuf) -> Result<faidx::Reader> {
    let reader = faidx::Reader::from_path(path).map_err(|e| e.to_string())?;
    Ok(reader)
}

fn is_bam_mapped(bam_header: &bam::Header) -> bool {
    // input is already sorted because it fails an index.
    // If it is mapped, the index needs the SQ tags to fetch data.
    for line in String::from_utf8(bam_header.to_bytes()).unwrap().lines() {
        if line.starts_with("@SQ") {
            return true;
        }
    }
    false
}

fn main() -> Result<()> {
    let params = get_cli_params();

    log::info!(
        "Running {}-{}",
        env!("CARGO_PKG_NAME"),
        *crate::cli::FULL_VERSION
    );
    let start_timer = time::Instant::now();

    let karyotype =
        Karyotype::new(&params.karyotype).unwrap_or_else(|err| handle_error_and_exit(err));

    let search_flank_len = params.flank_len;
    let output_flank_len = std::cmp::min(search_flank_len, 50);
    let sample_name = get_sample_name(&params.reads_path);

    let catalog_reader =
        open_catalog_reader(&params.repeats_path).unwrap_or_else(|err| handle_error_and_exit(err));
    let genome_reader = open_genome_reader(&params.genome_path)?;

    let all_loci = locus::get_loci(
        catalog_reader,
        &genome_reader,
        search_flank_len,
        &karyotype,
        params.genotyper,
    )
    .collect::<Result<Vec<_>>>()
    .unwrap_or_else(|err| handle_error_and_exit(err));

    let bam_header = get_bam_header(&params.reads_path)?;
    if !is_bam_mapped(&bam_header) {
        handle_error_and_exit("Input BAM is not mapped".into());
    }

    let mut vcf_writer = create_writer(&params.output_prefix, "vcf.gz", |path| {
        VcfWriter::new(path, &sample_name, &bam_header)
    })?;
    let mut bam_writer = create_writer(&params.output_prefix, "spanning.bam", |path| {
        BamWriter::new(path, bam_header)
    })?;

    log::info!("Starting job pool with {} threads...", params.num_threads);
    let pool: ThreadPool = ThreadPool::new(params.num_threads);
    let (sender, receiver) = channel();

    let writer_thread = thread::spawn(move || {
        for (locus, results) in &receiver {
            vcf_writer.write(&locus, &results);
            bam_writer.write(&locus, output_flank_len, &results);
        }
    });

    let reads_path = Arc::new(params.reads_path.clone());
    let workflow_params = Arc::new(workflows::Params {
        search_flank_len,
        min_read_qual: params.min_hifi_read_qual,
        max_depth: params.max_depth,
        aln_scoring: params.aln_scoring,
        min_flank_id_frac: params.min_flank_id_frac,
    });
    for locus in all_loci {
        let reads_path = reads_path.clone();
        let workflow_params = workflow_params.clone();
        let sender = sender.clone();

        pool.execute(move || {
            LOCAL.with(|local| {
                let mut bam = local.bam.borrow_mut();
                if bam.is_none() {
                    *bam = Some(
                        bam::IndexedReader::from_path(&*reads_path)
                            .expect("Failed to initialize bam file"),
                    );
                }
                let results = analyze_tr(&locus, &workflow_params, bam.as_mut().unwrap());
                match results {
                    Ok(results) => {
                        sender.send((locus, results)).unwrap();
                    }
                    Err(err) => {
                        eprintln!("Error occurred while analyzing: {}", err);
                    }
                }
            });
        });
    }
    pool.join();
    drop(sender);
    writer_thread.join().unwrap();
    log::info!("Total execution time: {:?}", start_timer.elapsed());
    log::info!("{} end", env!("CARGO_PKG_NAME"));
    Ok(())
}
