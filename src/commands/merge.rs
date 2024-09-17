use crate::cli::MergeArgs;
use crate::merge::vcf_processor::VcfProcessor;
use crate::utils::Result;
use std::time;

pub fn merge(args: MergeArgs) -> Result<()> {
    let start_timer = time::Instant::now();

    let mut vcf_processor = VcfProcessor::new(&args)?;

    if args.print_header {
        return Ok(());
    }

    vcf_processor.merge_variants()?;

    // TODO: If --output, --write-index is set and the output is compressed, index the file
    log::info!("Total execution time: {:.2?}", start_timer.elapsed());
    Ok(())
}
