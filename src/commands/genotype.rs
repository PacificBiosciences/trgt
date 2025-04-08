use crate::cli::GenotypeArgs;
use crate::trgt::{
    locus::{stream_loci_into_channel, Locus},
    workflows::{self, analyze_tr, LocusResult, Params},
    writers::{BamWriter, VcfWriter},
};
use crate::utils::{
    create_writer, get_bam_header, get_sample_name, is_bam_mapped, Karyotype, Result, TrgtScoring,
};
use crate::wfa_aligner::{
    AlignmentScope, Heuristic, MemoryModel, WFAligner, WFAlignerEdit, WFAlignerGapAffine,
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
    thread::{self},
};

#[derive(Debug, Clone)]
struct ThreadContextParams {
    flank_scoring: TrgtScoring,
    reads_path: PathBuf,
}

thread_local! {
    static CTX_PARAMS: RefCell<Option<ThreadContextParams>> = const { RefCell::new(None) };
}

fn create_thread_local_bam_reader() -> bam::IndexedReader {
    let path = CTX_PARAMS.with(|ctx_cell| {
        ctx_cell
            .borrow()
            .as_ref()
            .expect("Thread context parameters not initialized for BAM path")
            .reads_path
            .clone()
    });
    bam::IndexedReader::from_path(&path).unwrap_or_else(|e| {
        panic!(
            "Failed to initialize BAM reader for path {}: {}",
            path.display(),
            e
        )
    })
}

fn create_thread_local_ga_aligner_with_scoring() -> WFAligner {
    CTX_PARAMS.with(|ctx_cell| {
        let ctx = ctx_cell
            .borrow()
            .as_ref()
            .expect("Thread context parameters not initialized for WFA gap affine aligner")
            .clone();

        let scoring = &ctx.flank_scoring;
        let mut aligner = WFAlignerGapAffine::create_aligner_with_match(
            scoring.match_scr,
            scoring.mism_scr,
            scoring.gapo_scr,
            scoring.gape_scr,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        aligner.set_heuristic(Heuristic::None);
        aligner
    })
}

fn create_thread_local_ga_aligner() -> WFAligner {
    WFAlignerGapAffine::create_aligner_with_match(
        -1,
        1,
        5,
        1,
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
    )
}

fn create_thread_local_ed() -> WFAligner {
    WFAlignerEdit::create_aligner(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
}

thread_local! {
    // BAM Reader
    static THREAD_BAM_READER: RefCell<bam::IndexedReader> = RefCell::new(create_thread_local_bam_reader());
    // Flank locater aligner
    pub static THREAD_WFA_FLANK: RefCell<WFAligner> = RefCell::new(create_thread_local_ga_aligner_with_scoring());
    // Consensus aligner
    pub static THREAD_WFA_CONSENSUS: RefCell<WFAligner> = RefCell::new(create_thread_local_ga_aligner());
    // Edit distance
    pub static THREAD_WFA_ED: RefCell<WFAligner> = RefCell::new(create_thread_local_ed());
}

const CHANNEL_BUFFER_SIZE: usize = 2048;

pub fn trgt(args: GenotypeArgs) -> crate::utils::Result<()> {
    let karyotype = Karyotype::new(&args.karyotype)?;

    let bam_header = get_bam_header(&args.reads_path)?;
    if !is_bam_mapped(&bam_header) {
        return Err("Input BAM is not mapped".into());
    }

    let sample_name = args
        .sample_name
        .unwrap_or(get_sample_name(&args.reads_path, &bam_header)?);

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
        )
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
        min_flank_id_frac: args.min_flank_id_frac,
    });

    log::debug!(
        "Initializing thread pool with {} threads...",
        args.num_threads
    );

    let pool = initialize_thread_pool(
        args.num_threads,
        ThreadContextParams {
            flank_scoring: args.aln_scoring,
            reads_path: args.reads_path.clone(),
        },
    )?;
    pool.install(|| {
        receiver_locus
            .into_iter()
            .par_bridge()
            .for_each_with(&sender_result, |s, locus_result| match locus_result {
                Ok(locus) => process_locus(locus, &workflow_params, s),
                Err(err) => log::error!("Locus processing: {:#}", err),
            });
    });

    // Clean-up
    drop(sender_result);
    writer_thread.join().expect("Writer thread panicked");
    log::trace!("Writer thread finished");
    match locus_stream_thread
        .join()
        .expect("Locus stream thread panicked")
    {
        Ok(_) => log::trace!("Locus stream thread finished"),
        Err(e) => log::error!("Locus streaming failed: {}", e),
    }

    Ok(())
}

fn process_locus(
    locus: Locus,
    workflow_params: &Arc<Params>,
    sender_result: &Sender<(Locus, LocusResult)>,
) {
    THREAD_BAM_READER.with(|reader_cell| {
        let mut reader = reader_cell.borrow_mut();
        match analyze_tr(&locus, workflow_params, &mut reader) {
            Ok(results) => {
                if let Err(e) = sender_result.send((locus, results)) {
                    log::error!("Failed to send locus result to writer thread: {}", e);
                }
            }
            Err(err) => {
                log::error!("Error analyzing locus {}: {}", locus.id, err);
            }
        }
    });
}

fn initialize_thread_pool(
    num_threads: usize,
    thread_context: ThreadContextParams,
) -> Result<rayon::ThreadPool> {
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .thread_name(|i| format!("trgt-{}", i))
        .start_handler(move |_thread_index| {
            CTX_PARAMS.with(|cell| {
                *cell.borrow_mut() = Some(thread_context.clone());
            });
            log::trace!("Initialized thread {:?}", std::thread::current().id());
        })
        .exit_handler(|_thread_index| {
            CTX_PARAMS.with(|cell| {
                *cell.borrow_mut() = None;
            });
        })
        .build()
        .map_err(|e| format!("Failed to initialize thread pool: {}", e))
}
