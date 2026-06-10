mod bam;
mod catalog;
mod genome;
mod vcf;

pub use bam::{
    create_bam_reader_with_options, open_bam_reader, open_genotyper_bam_reader,
    open_genotyper_bam_reader_with_thread_pool, BamOpenOptions, ReferencePolicy,
};
pub use catalog::open_catalog_reader;
pub use genome::{check_missing_faidx, open_genome_reader};
pub use vcf::open_vcf_reader;
