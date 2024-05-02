use crate::cli::ValidateArgs;
use crate::trgt::locus::get_loci;
use crate::utils::{open_catalog_reader, open_genome_reader, Genotyper, Karyotype, Result};

pub fn validate(args: ValidateArgs) -> Result<()> {
    let catalog_reader = open_catalog_reader(&args.repeats_path)?;
    let genome_reader = open_genome_reader(&args.genome_path)?;
    let mut error_count = 0;
    let mut success_count = 0;
    let mut motifs_lengths = Vec::new();
    let mut tr_lengths = Vec::new();

    for result in get_loci(
        catalog_reader,
        genome_reader,
        Karyotype::new("XY").unwrap(),
        args.flank_len,
        Genotyper::Size,
    ) {
        match result {
            Ok(locus) => {
                motifs_lengths.push(locus.motifs.len());
                tr_lengths.push(locus.tr.len());
                success_count += 1
            }
            Err(e) => {
                log::error!("{}", e);
                error_count += 1;
            }
        }
    }

    let motifs_stats = calculate_stats(&motifs_lengths);
    let tr_stats = calculate_stats(&tr_lengths);

    let total = success_count + error_count;
    let success_percentage = (success_count as f64 / total as f64) * 100.0;
    let error_percentage = (error_count as f64 / total as f64) * 100.0;

    log::info!(
        "Motifs per Locus - Range: [{},{}], Median: {:.2}, Mean: {:.2}, StdDev: {:.2}",
        motifs_stats.min,
        motifs_stats.max,
        motifs_stats.mean,
        motifs_stats.median,
        motifs_stats.std_dev
    );
    log::info!(
        "TR Lengths - Range: [{},{}], Median: {:.2}, Mean: {:.2}, StdDev: {:.2}",
        tr_stats.min,
        tr_stats.max,
        tr_stats.mean,
        tr_stats.median,
        tr_stats.std_dev
    );

    match error_count {
        0 => log::info!("Validation successful. Loci pass={}", success_count),
        _ => log::info!(
            "Validation failed. Loci pass={} ({:.2}%), fail={} ({:.2}%)",
            success_count,
            success_percentage,
            error_count,
            error_percentage
        ),
    }

    Ok(())
}

fn calculate_stats(data: &[usize]) -> Stats {
    let mut sorted = data.to_vec();
    sorted.sort_unstable();
    let len = sorted.len();
    let median = if len % 2 == 0 {
        (sorted[len / 2 - 1] + sorted[len / 2]) as f64 / 2.0
    } else {
        sorted[len / 2] as f64
    };
    let sum: usize = sorted.iter().sum();
    let mean = sum as f64 / len as f64;
    let std_dev = (sorted
        .iter()
        .map(|&x| (x as f64 - mean).powi(2))
        .sum::<f64>()
        / len as f64)
        .sqrt();
    Stats {
        min: *sorted.first().unwrap_or(&0),
        max: *sorted.last().unwrap_or(&0),
        mean,
        median,
        std_dev,
    }
}

struct Stats {
    min: usize,
    max: usize,
    mean: f64,
    median: f64,
    std_dev: f64,
}
