use crate::cli::MergeArgs;
use crate::merge::vcf_processor::VcfProcessor;
use crate::utils::Result;

pub fn merge(args: MergeArgs) -> Result<()> {
    let vcf_paths = args.process_vcf_paths()?;
    let mut vcf_processor = VcfProcessor::new(&args, vcf_paths)?;

    if args.print_header {
        return Ok(());
    }

    vcf_processor.merge_variants()?;

    // TODO: If --output, --write-index is set and the output is compressed, index the file
    Ok(())
}
