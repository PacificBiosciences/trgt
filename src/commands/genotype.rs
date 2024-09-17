use crate::cli::GenotypeArgs;
use crate::trgt::{
    locus::{stream_loci_into_channel, Locus},
    workflows::{self, analyze_tr, LocusResult, Params},
    writers::{BamWriter, VcfWriter},
};
use crate::utils::{
    create_writer, get_bam_header, get_sample_name, is_bam_mapped, open_catalog_reader,
    open_genome_reader, Karyotype, Result,
};
use crossbeam_channel::{bounded, Sender};
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    ThreadPoolBuilder,
};
use rust_htslib::bam;
use std::{
    cell::RefCell,
    path::PathBuf,
    sync::Arc,
    thread::{self, available_parallelism},
    time,
};

struct ThreadLocalData {
    bam: RefCell<Option<bam::IndexedReader>>,
}

thread_local! {
    static LOCAL_BAM_READER: ThreadLocalData = const { ThreadLocalData {
        bam: RefCell::new(None),
    } };
}

const CHANNEL_BUFFER_SIZE: usize = 2048;

pub fn trgt(args: GenotypeArgs) -> Result<()> {
    let start_timer = time::Instant::now();

    let karyotype = Karyotype::new(&args.karyotype)?;

    let bam_header = get_bam_header(&args.reads_path)?;
    if !is_bam_mapped(&bam_header) {
        return Err("Input BAM is not mapped".into());
    }

    let sample_name = args
        .sample_name
        .unwrap_or(get_sample_name(&args.reads_path, &bam_header)?);

    // TODO: Find a nicer solution, this is still necessary since we do want to check if we can open the catalog and genome reader successfully; note that we also open these inside stream_loci_into_channel...
    open_catalog_reader(&args.repeats_path)?;
    open_genome_reader(&args.genome_path)?;

    let mut vcf_writer = create_writer(&args.output_prefix, "vcf.gz", |path| {
        VcfWriter::new(path, &sample_name, &bam_header)
    })?;

    let output_flank_len = std::cmp::min(args.flank_len, args.output_flank_len);
    let bam_writer = if !args.disable_bam_output {
        Some(create_writer(
            &args.output_prefix,
            "spanning.bam",
            |path| BamWriter::new(path, bam_header, output_flank_len),
        )?)
    } else {
        None
    };

    let (sender_locus, receiver_locus) = bounded(CHANNEL_BUFFER_SIZE);
    let locus_stream_thread = thread::spawn(move || {
        stream_loci_into_channel(
            &args.repeats_path,
            &args.genome_path,
            args.flank_len,
            args.genotyper,
            &karyotype,
            sender_locus,
        );
    });

    let (sender_result, receiver_result) = bounded(CHANNEL_BUFFER_SIZE);
    let writer_thread = thread::spawn(move || {
        if let Some(mut bam_writer) = bam_writer {
            for (locus, results) in &receiver_result {
                vcf_writer.write(&locus, &results);
                bam_writer.write(&locus, &results);
            }
        } else {
            for (locus, results) in &receiver_result {
                vcf_writer.write(&locus, &results);
            }
        }
    });

    let workflow_params = Arc::new(workflows::Params {
        search_flank_len: args.flank_len,
        min_read_qual: args.min_hifi_read_qual,
        max_depth: args.max_depth,
        aln_scoring: args.aln_scoring,
        min_flank_id_frac: args.min_flank_id_frac,
    });

    if args.num_threads == 1 {
        log::debug!("Single-threaded mode");
        for locus in receiver_locus {
            match locus {
                Ok(locus) => {
                    process_locus(locus, &args.reads_path, &workflow_params, &sender_result)
                }
                Err(err) => log::error!("Locus Processing: {:#}", err),
            }
        }
    } else {
        log::debug!(
            "Multi-threaded mode: estimated available cores: {}",
            available_parallelism().unwrap().get()
        );
        log::info!("Starting job pool with {} threads...", args.num_threads);
        let pool = initialize_thread_pool(args.num_threads)?;

        pool.install(|| {
            receiver_locus
                .into_iter()
                .par_bridge()
                .for_each_with(&sender_result, |s, locus| match locus {
                    Ok(locus) => process_locus(locus, &args.reads_path, &workflow_params, s),
                    Err(err) => log::error!("Locus Processing: {:#}", err),
                });
        });
    }
    drop(sender_result);
    writer_thread.join().unwrap();
    locus_stream_thread.join().unwrap();
    log::info!("Total execution time: {:.2?}", start_timer.elapsed());
    Ok(())
}

fn process_locus(
    locus: Locus,
    reads_path: &PathBuf,
    workflow_params: &Arc<Params>,
    sender_result: &Sender<(Locus, LocusResult)>,
) {
    LOCAL_BAM_READER.with(|local| {
        let mut bam = local.bam.borrow_mut();
        if bam.is_none() {
            *bam = Some(
                bam::IndexedReader::from_path(reads_path).expect("Failed to initialize bam file"),
            );
        }
        match analyze_tr(&locus, workflow_params, bam.as_mut().unwrap()) {
            Ok(results) => {
                sender_result.send((locus, results)).unwrap();
            }
            Err(err) => {
                log::error!("Error occurred while analyzing: {}", err);
            }
        }
    });
}

fn initialize_thread_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| format!("Failed to initialize thread pool: {}", e))
}
