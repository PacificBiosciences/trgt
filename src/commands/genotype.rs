use crate::cli::GenotypeArgs;
use crate::trgt::{
    genotype_reader::{BamRead, BamSourceReader, MultiBamReader, SingleBamReader},
    locus::{stream_locus_groups_into_channel, Locus},
    locus_group::{extract_reads_for_group, LocusGroup},
    workflows::{self, analyze_tr, LocusResult, Params},
    writers::{BamWriter, VcfWriter},
};
use crate::utils::{
    create_writer,
    readers::{open_genotyper_bam_reader, open_genotyper_bam_reader_with_thread_pool},
    validate_bam_sources, BamSources, InputSource, Result, TrgtScoring,
};
use crate::wfaligner::{AlignmentScope, Heuristic, MemoryModel, WFAligner};
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::{
    iter::{ParallelBridge, ParallelIterator},
    ThreadPoolBuilder,
};
use rust_htslib::tpool::ThreadPool;
use std::{cell::RefCell, sync::Arc, thread};

#[derive(Debug, Clone)]
struct ThreadContextParams {
    flank_scoring: TrgtScoring,
}

thread_local! {
    static THREAD_CONTEXT: RefCell<Option<ThreadContextParams>> = const { RefCell::new(None) };
}

fn create_thread_local_ga_aligner_with_scoring() -> WFAligner {
    THREAD_CONTEXT.with(|context_cell| {
        let context = context_cell
            .borrow()
            .as_ref()
            .expect("Thread context parameters not initialized for WFA gap affine aligner")
            .clone();
        let scoring = &context.flank_scoring;

        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(scoring.mism_scr, scoring.gapo_scr, scoring.gape_scr)
            .with_heuristic(Heuristic::None)
            .build()
    })
}

fn create_thread_local_ga_aligner() -> WFAligner {
    WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(2, 5, 1)
        .build()
}

fn create_thread_local_ed() -> WFAligner {
    WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
        .edit()
        .build()
}

thread_local! {
    // Flank locater aligner
    pub static THREAD_WFA_FLANK: RefCell<WFAligner> = RefCell::new(create_thread_local_ga_aligner_with_scoring());
    // Consensus aligner
    pub static THREAD_WFA_CONSENSUS: RefCell<WFAligner> = RefCell::new(create_thread_local_ga_aligner());
    // Edit distance
    pub static THREAD_WFA_ED: RefCell<WFAligner> = RefCell::new(create_thread_local_ed());
}

pub fn trgt(args: GenotypeArgs) -> Result<()> {
    let thread_context = ThreadContextParams {
        flank_scoring: args.aln_scoring,
    };

    let bam_sources = args.read_sources();

    let (bam_header, sample_name) =
        validate_bam_sources(&bam_sources, args.sample_name.as_deref())?;

    // TODO: A thread pool for the writers can be shared between the VCF and BAM writers
    let output_type = args.output_type.unwrap_or_default();
    let output_suffix = output_type.output_suffix();
    let mut vcf_writer = create_writer(&args.output_prefix, output_suffix, |path| {
        VcfWriter::new(path, &sample_name, &bam_header, Some(output_type))
    })?;

    let output_flank_len = std::cmp::min(args.flank_len, args.output_flank_len);
    let bam_writer = if !args.disable_bam_output {
        let n_writer_threads = args.bam_compression_threads.unwrap_or({
            match args.num_threads {
                n if n < 12 => 1,
                _ => 2,
            }
        });
        Some(create_writer(
            &args.output_prefix,
            "spanning.bam",
            |path| BamWriter::new(path, bam_header, output_flank_len, n_writer_threads),
        )?)
    } else {
        None
    };

    let workflow_params = Arc::new(workflows::Params {
        search_flank_len: args.flank_len,
        min_read_qual: args.min_hifi_read_qual,
        max_depth: args.max_depth,
        min_flank_id_frac: args.min_flank_id_frac,
        skip_phase_annotation: args.skip_phase_annotation,
    });

    let num_fetcher_threads = args
        .fetcher_threads
        .unwrap_or_else(|| (args.num_threads / 2).clamp(1, 8));

    let group_channel_buffer_size = args
        .group_channel_buffer_size
        .unwrap_or(num_fetcher_threads);

    with_analysis_thread_pool(args.num_threads, thread_context, |pool| {
        thread::scope(|scope| -> Result<()> {
            // Stage 1 -> 2: Locus group
            let (locus_group_tx, locus_group_rx) = bounded(group_channel_buffer_size);
            // Stage 2 -> 3: Loci with reads
            let (populated_locus_tx, populated_locus_rx) = bounded(args.locus_channel_buffer_size);
            // Stage 3 -> 4: Final results for writer
            let (analysis_result_tx, analysis_result_rx) = bounded(args.result_channel_buffer_size);

            // Stage 4: Spawn the writer thread to write BAM/VCF output (IO: Write BAM and VCF)
            let writer_thread = spawn_scoped_stage(scope, "writer", move || {
                if let Some(mut bam_writer) = bam_writer {
                    for (locus, results) in &analysis_result_rx {
                        vcf_writer.write(&locus, &results);
                        bam_writer.write(&locus, &results);
                    }
                } else {
                    for (locus, results) in &analysis_result_rx {
                        vcf_writer.write(&locus, &results);
                    }
                }
                Ok(())
            })?;

            // Stage 2: Populate LocusGroup with reads (IO: Read BAM file)
            let fetcher_threads = spawn_fetcher_threads(
                scope,
                &bam_sources,
                &args.genome_src,
                num_fetcher_threads,
                args.decompression_threads,
                &workflow_params,
                locus_group_rx,
                &populated_locus_tx,
            )?;
            drop(populated_locus_tx);

            // Stage 1: Grouping loci into LocusGroup's (IO: Read repeat catalog and reference genome)
            let locus_grouper_thread = {
                // Move the only group sender into the grouper so its exit closes stage 1
                spawn_scoped_stage(scope, "locus-grouper", move || {
                    stream_locus_groups_into_channel(
                        &args.repeats_src,
                        &args.genome_src,
                        args.flank_len,
                        args.genotyper,
                        &args.karyotype,
                        locus_group_tx,
                        args.max_group_span,
                        args.max_group_size,
                    )
                    .map_err(|e| format!("Locus grouping stage failed: {e}"))
                })?
            };

            // Stage 3: Worker pool for genotyping TRs
            let analysis_sender_for_workers = analysis_result_tx.clone();
            let analysis_result = pool.install(|| {
                populated_locus_rx
                    .into_iter()
                    .par_bridge()
                    .try_for_each_with(analysis_sender_for_workers, |s, locus_result| {
                        match locus_result {
                            Ok(locus) => process_locus(locus, &workflow_params, s),
                            Err(err) => {
                                log::error!("Locus processing: {err:#}");
                                Ok(())
                            }
                        }
                    })
            });
            // Worker sender clones are gone after analysis, close the writer input
            drop(analysis_result_tx);

            let writer_result = writer_thread
                .join()
                .unwrap_or_else(|_| Err("Writer thread panicked.".into()));
            let mut stage_result = match writer_result {
                Ok(()) => {
                    log::trace!("Writer thread has finished.");
                    analysis_result
                }
                Err(writer_error) => {
                    if let Err(analysis_error) = analysis_result {
                        log::error!(
                            "Stopped analysis stage after writer failure: {analysis_error}"
                        );
                    }
                    Err(writer_error)
                }
            };

            let locus_grouper_result = locus_grouper_thread
                .join()
                .unwrap_or_else(|_| Err("Locus grouper thread panicked.".into()));
            match locus_grouper_result {
                Ok(()) => log::trace!("Locus grouper thread has finished."),
                Err(e) if stage_result.is_ok() => stage_result = Err(e),
                Err(_) => {}
            }

            for handle in fetcher_threads {
                let fetcher_result = handle
                    .join()
                    .unwrap_or_else(|_| Err("Fetcher thread panicked.".into()));
                match fetcher_result {
                    Ok(()) => {}
                    Err(e) if stage_result.is_ok() => stage_result = Err(e),
                    Err(_) => {}
                }
            }
            stage_result?;
            log::trace!("All fetcher threads have finished.");
            Ok(())
        })
    })?;

    log::debug!("Analysis complete. All pipeline threads shut down successfully.");
    Ok(())
}

fn run_fetcher_thread<R: BamRead>(
    mut bam_reader: R,
    locus_group_rx: Receiver<Result<LocusGroup>>,
    populated_locus_tx: Sender<Result<Locus>>,
    workflow_params: Arc<Params>,
) {
    for group_result in locus_group_rx {
        let group = match group_result {
            Ok(group) => group,
            Err(e) => {
                if populated_locus_tx.send(Err(e)).is_err() {
                    return;
                }
                continue;
            }
        };

        match extract_reads_for_group(&mut bam_reader, group, &workflow_params) {
            Ok(populated_loci) => {
                for locus in populated_loci {
                    if populated_locus_tx.send(Ok(locus)).is_err() {
                        return;
                    }
                }
            }
            Err(e) => {
                if populated_locus_tx.send(Err(e)).is_err() {
                    return;
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn spawn_fetcher_threads<'scope>(
    scope: &'scope thread::Scope<'scope, '_>,
    bam_sources: &BamSources,
    genome_src: &InputSource,
    num_fetcher_threads: usize,
    decompression_threads: usize,
    workflow_params: &Arc<Params>,
    locus_group_rx: Receiver<Result<LocusGroup>>,
    populated_locus_tx: &Sender<Result<Locus>>,
) -> Result<Vec<thread::ScopedJoinHandle<'scope, Result<()>>>> {
    let mut threads = Vec::with_capacity(num_fetcher_threads);
    let read_source_count = bam_sources.as_slice().len();
    let open_reader_count = num_fetcher_threads * read_source_count;
    log::debug!(
        "Initializing read fetcher pool with {num_fetcher_threads} fetcher thread(s), {read_source_count} read source(s), implied open-reader count: {open_reader_count}, and decompression pool size per fetcher: {decompression_threads} thread(s)."
    );
    match bam_sources.as_slice() {
        [] => return Err("At least one BAM source is required".to_string()),
        [source] => {
            for i in 0..num_fetcher_threads {
                let reads_src = source.clone();
                let genome_src = genome_src.clone();
                let locus_group_rx = locus_group_rx.clone();
                let populated_locus_tx = populated_locus_tx.clone();
                let workflow_params = workflow_params.clone();
                let handle =
                    spawn_scoped_stage(scope, &format!("fetcher-{i}"), move || -> Result<()> {
                        let bam_reader = open_single_fetcher_bam_reader(
                            &reads_src,
                            &genome_src,
                            decompression_threads,
                        )
                        .map_err(|e| {
                            format!("Thread fetcher-{i} failed to open BAM reader: {e}")
                        })?;
                        run_fetcher_thread(
                            bam_reader,
                            locus_group_rx,
                            populated_locus_tx,
                            workflow_params,
                        );
                        Ok(())
                    })?;

                threads.push(handle);
            }
        }
        _ => {
            for i in 0..num_fetcher_threads {
                let bam_sources = bam_sources.clone();
                let genome_src = genome_src.clone();
                let locus_group_rx = locus_group_rx.clone();
                let populated_locus_tx = populated_locus_tx.clone();
                let workflow_params = workflow_params.clone();
                let handle =
                    spawn_scoped_stage(scope, &format!("fetcher-{i}"), move || -> Result<()> {
                        let bam_reader = open_multi_fetcher_bam_reader(
                            &bam_sources,
                            &genome_src,
                            decompression_threads,
                        )
                        .map_err(|e| {
                            format!("Thread fetcher-{i} failed to open BAM reader: {e}")
                        })?;
                        run_fetcher_thread(
                            bam_reader,
                            locus_group_rx,
                            populated_locus_tx,
                            workflow_params,
                        );
                        Ok(())
                    })?;

                threads.push(handle);
            }
        }
    }
    Ok(threads)
}

fn open_single_fetcher_bam_reader(
    reads_src: &InputSource,
    genome_src: &InputSource,
    decompression_threads: usize,
) -> Result<SingleBamReader> {
    open_genotyper_bam_reader(reads_src, genome_src, decompression_threads)
        .map(SingleBamReader::new)
}

fn open_multi_fetcher_bam_reader(
    bam_sources: &BamSources,
    genome_src: &InputSource,
    decompression_threads: usize,
) -> Result<MultiBamReader> {
    let thread_pool = {
        let thread_count = u32::try_from(decompression_threads).map_err(|_| {
            format!("Decompression thread count {decompression_threads} exceeds u32::MAX")
        })?;
        ThreadPool::new(thread_count)
            .map_err(|e| format!("Failed to initialize BAM decompression thread pool: {e}"))?
    };
    let source_readers = bam_sources
        .iter_labeled()
        .map(|(source, label)| {
            open_genotyper_bam_reader_with_thread_pool(source, genome_src, &thread_pool)
                .map(|reader| BamSourceReader::new(label, reader))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(MultiBamReader::with_tpool(source_readers, thread_pool))
}

fn process_locus(
    locus: Locus,
    workflow_params: &Params,
    sender_result: &Sender<(Locus, LocusResult)>,
) -> Result<()> {
    let (locus, result) = analyze_tr(locus, workflow_params);
    match result {
        Ok(results) => sender_result
            .send((locus, results))
            .map_err(|e| format!("Failed to send locus result to writer thread: {e}")),
        Err(err) => {
            log::error!(
                "Error analyzing locus {locus_id}: {err}",
                locus_id = locus.id
            );
            Ok(())
        }
    }
}

fn spawn_scoped_stage<'scope, 'env, F>(
    scope: &'scope thread::Scope<'scope, 'env>,
    thread_name: &str,
    f: F,
) -> Result<thread::ScopedJoinHandle<'scope, Result<()>>>
where
    F: FnOnce() -> Result<()> + Send + 'scope,
{
    thread::Builder::new()
        .name(thread_name.to_owned())
        .spawn_scoped(scope, move || {
            let thread_current = thread::current();
            let thread_name = thread_current.name().unwrap_or("<unnamed>");
            log::trace!("Starting thread {thread_name}.");
            let result = f();
            match &result {
                Ok(()) => {
                    log::trace!("Exiting thread: {thread_name} finished.");
                }
                Err(e) => {
                    log::error!("Exiting thread: {thread_name} failed: {e}");
                }
            }
            result
        })
        .map_err(|e| format!("Failed to spawn {thread_name}: {e}"))
}

fn with_analysis_thread_pool<F, R>(
    num_threads: usize,
    thread_context: ThreadContextParams,
    with_pool: F,
) -> Result<R>
where
    F: FnOnce(&rayon::ThreadPool) -> Result<R>,
{
    log::debug!("Initializing analysis thread pool with {num_threads} threads.");
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .thread_name(|i| format!("trgt-analysis-{i}"))
        .start_handler(move |_thread_index| {
            let thread_current = std::thread::current();
            let thread_name = thread_current.name().unwrap_or("<unnamed>");
            log::trace!("Initializing thread: {thread_name}.");

            THREAD_CONTEXT.with(|cell| {
                *cell.borrow_mut() = Some(thread_context.clone());
            });
            THREAD_WFA_FLANK.with(|_| {});
            THREAD_WFA_CONSENSUS.with(|_| {});
            THREAD_WFA_ED.with(|_| {});
        })
        .exit_handler(|_thread_index| {
            let thread_current = std::thread::current();
            let thread_name = thread_current.name().unwrap_or("<unnamed>");

            THREAD_CONTEXT.with(|cell| {
                *cell.borrow_mut() = None;
            });
            log::trace!("Exiting thread: {thread_name}.");
        })
        .build_scoped(|thread| thread.run(), with_pool)
        .map_err(|e| format!("Failed to initialize thread pool: {e}"))?
}
