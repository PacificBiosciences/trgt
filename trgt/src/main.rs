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

use cli::{get_cli_params, handle_error_and_exit};
use flate2::read::GzDecoder;
use karyotype::Karyotype;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::cell::RefCell;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, Read as ioRead};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::{thread, time};
use threadpool::ThreadPool;
use workflows::analyze_tr;
use writers::{BamWriter, VcfWriter};
mod cli;
mod cluster;
mod faidx;
mod genotype;
mod hmm;
mod karyotype;
mod locate;
mod locus;
mod reads;
mod utils;
mod workflows;
mod writers;

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
    let bam = bam::IndexedReader::from_path(bam_path)
        .map_err(|e| format!("Failed to create bam reader: {}", e))?;
    Ok(bam::Header::from_template(bam.header()))
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

fn get_sample_name(reads_path: &PathBuf) -> Result<String> {
    let bam_header = get_bam_header(reads_path)?;

    let header_hashmap = bam_header.to_hashmap();
    let mut sample_names = HashSet::new();

    if let Some(rg_fields) = header_hashmap.get("RG") {
        for rg_field in rg_fields {
            if let Some(sample_name) = rg_field.get("SM") {
                sample_names.insert(sample_name.to_owned());
            }
        }
    }

    match sample_names.len() {
        1 => return Ok(sample_names.into_iter().next().unwrap()),
        0 => log::warn!("No sample names found"),
        _ => log::warn!("Multiple sample names found"),
    };

    let sample = reads_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("Invalid reads file name")?
        .to_string();

    Ok(sample)
}

fn create_writer<T, F>(output_prefix: &str, output_suffix: &str, f: F) -> Result<T>
where
    F: FnOnce(&str) -> Result<T>,
{
    let output_path = format!("{}.{}", output_prefix, output_suffix);
    f(&output_path)
}

fn open_catalog_reader(path: &PathBuf) -> Result<BufReader<Box<dyn ioRead>>> {
    fn is_gzipped(path: &Path) -> bool {
        let path_str = path.to_string_lossy();
        let formats = [".gz", ".gzip", ".GZ", ".GZIP"];
        formats.iter().any(|format| path_str.ends_with(*format))
    }
    let file = File::open(path).map_err(|e| e.to_string())?;
    if is_gzipped(path) {
        let gz_decoder = GzDecoder::new(file);
        if gz_decoder.header().is_some() {
            Ok(BufReader::new(Box::new(gz_decoder)))
        } else {
            Err(format!("Invalid gzip header: {}", path.to_string_lossy()))
        }
    } else {
        Ok(BufReader::new(Box::new(file)))
    }
}

fn open_genome_reader(path: &Path) -> Result<faidx::Reader> {
    let extension = path.extension().unwrap().to_str().unwrap();
    let fai_path = path.with_extension(extension.to_owned() + ".fai");
    if !fai_path.exists() {
        return Err(format!(
            "Reference index file not found: {}. Create it using 'samtools faidx {}'",
            fai_path.display(),
            path.display()
        ));
    }
    faidx::Reader::from_path(path).map_err(|e| e.to_string())
}

fn main() {
    if let Err(e) = run_trgt() {
        handle_error_and_exit(e);
    }
}

fn run_trgt() -> Result<()> {
    let params = get_cli_params();

    log::info!(
        "Running {}-{}",
        env!("CARGO_PKG_NAME"),
        *crate::cli::FULL_VERSION
    );
    let start_timer = time::Instant::now();

    let karyotype = Karyotype::new(&params.karyotype)?;

    let sample_name = params
        .sample_name
        .unwrap_or(get_sample_name(&params.reads_path)?);

    let catalog_reader = open_catalog_reader(&params.repeats_path)?;
    let genome_reader = open_genome_reader(&params.genome_path)?;

    let all_loci = locus::get_loci(
        catalog_reader,
        &genome_reader,
        karyotype,
        params.flank_len,
        params.genotyper,
    )
    .collect::<Result<Vec<_>>>()?;

    let bam_header = get_bam_header(&params.reads_path)?;
    if !is_bam_mapped(&bam_header) {
        handle_error_and_exit("Input BAM is not mapped".into());
    }

    let mut vcf_writer = create_writer(&params.output_prefix, "vcf.gz", |path| {
        VcfWriter::new(path, &sample_name, &bam_header)
    })?;

    let output_flank_len = std::cmp::min(params.flank_len, params.output_flank_len);
    let mut bam_writer = create_writer(&params.output_prefix, "spanning.bam", |path| {
        BamWriter::new(path, bam_header, output_flank_len)
    })?;

    log::info!("Starting job pool with {} threads...", params.num_threads);
    let pool: ThreadPool = ThreadPool::new(params.num_threads);
    let (sender, receiver) = channel();

    let writer_thread = thread::spawn(move || {
        for (locus, results) in &receiver {
            vcf_writer.write(&locus, &results);
            bam_writer.write(&locus, &results);
        }
    });

    let reads_path = Arc::new(params.reads_path.clone());
    let workflow_params = Arc::new(workflows::Params {
        search_flank_len: params.flank_len,
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
                        log::error!("Error occurred while analyzing: {}", err);
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
