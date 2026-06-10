use crate::utils::{readers::check_missing_faidx, InputSource, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::tpool::ThreadPool;
use std::path::Path;

pub enum ReferencePolicy {
    Never,
    Always,
    AutoForCram,
}

pub struct BamOpenOptions<'a> {
    pub threads: Option<usize>,
    pub thread_pool: Option<&'a ThreadPool>,
    pub reference: Option<&'a InputSource>,
    pub ref_policy: ReferencePolicy,
}

impl Default for BamOpenOptions<'_> {
    fn default() -> Self {
        Self {
            threads: None,
            thread_pool: None,
            reference: None,
            ref_policy: ReferencePolicy::Never,
        }
    }
}

fn is_cram_source(src: &InputSource) -> bool {
    match src {
        InputSource::Local(p) => has_ext(p, "cram"),
        InputSource::Remote(r) => r.url().path().to_ascii_lowercase().ends_with(".cram"),
    }
}

fn has_ext(p: &Path, ext: &str) -> bool {
    p.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case(ext))
        .unwrap_or(false)
}

pub fn create_bam_reader_with_options(
    reads_src: &InputSource,
    opts: &BamOpenOptions<'_>,
) -> Result<bam::IndexedReader> {
    reads_src.preflight_checks()?;
    let mut reader = match reads_src {
        InputSource::Local(p) => bam::IndexedReader::from_path(p),
        InputSource::Remote(r) => bam::IndexedReader::from_url(r.url()),
    }
    .map_err(|e| reads_src.format_error("Failed to create BAM reader from", e))?;

    match (opts.threads, opts.thread_pool) {
        (Some(_), Some(_)) => {
            return Err(
                "Cannot set both decompression threads and a decompression thread pool".to_string(),
            );
        }
        (Some(n), None) => {
            if let Err(e) = reader.set_threads(n) {
                log::warn!("Failed to set decompression threads: {e}");
            }
        }
        (None, Some(thread_pool)) => {
            reader
                .set_thread_pool(thread_pool)
                .map_err(|e| format!("Failed to set decompression thread pool: {e}"))?;
        }
        (None, None) => {}
    }

    let need_ref = match opts.ref_policy {
        ReferencePolicy::Never => false,
        ReferencePolicy::Always => true,
        ReferencePolicy::AutoForCram => is_cram_source(reads_src),
    };
    if need_ref {
        if let Some(ref_src) = opts.reference {
            if let InputSource::Local(p) = ref_src {
                check_missing_faidx(p)?;
                if let Err(e) = reader.set_reference(p) {
                    log::warn!("Failed to set reference {}: {e}", p.display());
                }
            } else {
                log::warn!("Reference is remote; skipping set_reference()");
            }
        }
    }
    Ok(reader)
}

pub fn open_genotyper_bam_reader(
    reads_src: &InputSource,
    genome_src: &InputSource,
    decompression_threads: usize,
) -> Result<bam::IndexedReader> {
    let opts = BamOpenOptions {
        threads: Some(decompression_threads),
        thread_pool: None,
        reference: Some(genome_src),
        ref_policy: ReferencePolicy::AutoForCram,
    };
    create_bam_reader_with_options(reads_src, &opts)
}

pub fn open_genotyper_bam_reader_with_thread_pool(
    reads_src: &InputSource,
    genome_src: &InputSource,
    thread_pool: &ThreadPool,
) -> Result<bam::IndexedReader> {
    let opts = BamOpenOptions {
        threads: None,
        thread_pool: Some(thread_pool),
        reference: Some(genome_src),
        ref_policy: ReferencePolicy::AutoForCram,
    };
    create_bam_reader_with_options(reads_src, &opts)
}

pub fn open_bam_reader(reads_src: &InputSource, threads: usize) -> Result<bam::IndexedReader> {
    let opts = BamOpenOptions {
        threads: Some(threads),
        ..Default::default()
    };
    create_bam_reader_with_options(reads_src, &opts)
}
