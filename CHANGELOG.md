
# Changelog

Any changes to TRGT are noted below:

## [5.1.0]

- **Important changes**:
  - Added the `--fail-reads` parameter to `trgt genotype`. It is now possible to add a fail reads BAM as a separate aligned BAM in addition to the primary `--reads` BAM. This is useful for PureTarget and other workflows where fail and HiFi reads are aligned separately, since they no longer need to be merged before running TRGT. Adding `--fail-reads` also disables the default `rq` filter unless `--min-read-quality` is set explicitly.
- Added BCF output support to `trgt genotype` through `-O, --output-type` (`u`/`b` for uncompressed/compressed BCF, `v`/`z` for uncompressed/compressed VCF).

## [5.0.0]

- **Important changes**:
  - Extended support for phased BAM inputs (from tool such as [HiPhase](https://github.com/PacificBiosciences/HiPhase)). TRGT can now output phased genotypes (e.g. `1|0`) and phase sets using HP/PS tags, when sufficient haplotype support is available.
    - **Note:** the FORMAT field order remains unchanged after phasing, any downstream analysis that assumes positional indexing on genotypes must account for phasing.
  - More robust and optimized VCF merging (`merge` subcommand).
    - Added a streaming merge mode which avoids loading tabix/CSI indices. This mode is especially useful for merging very large cohorts and can be enabled with the `--no-index` flag. It significantly reduces memory requirements but it requires that all inputs have the exact same contig ordering.
    - Compressed VCF/BCF output can now be automatically indexed after merging by setting the `--write-index` flag.
    - Added the `--threads` flag to set the number of (de)compression threads for input/output VCF files.

## [4.1.0]

- **Important changes**:
  - Experimental remote file support for most TRGT inputs (BAM/CRAM, reference FASTA, repeat catalogs, and VCF/BCF). See [remote file documentation](docs/remote_files.md) for more details.
  - Improved consensus generation for long repeat expansions, as part of an ongoing refinement effort.
- CLI: improved terminal color detection, and added a `--color` flag.

## [4.0.0]

- **Important changes**:
  - Implemented read streaming to reduce I/O operations, achieving 2-4x speed improvements for WGS data (speed-up varies by catalog density). The input catalog must be sorted by contig and start to achieve high performance.
  - Added `trgt deepdive` subcommand for locus-specific downstream analysis. This command realigns reads to their consensus sequences and outputs FASTA (consensus sequences), BAM (realigned reads), and BED (annotation) files for detailed repeat analysis in tools like IGV. See [deepdive documentation](docs/deepdive.md) for more details.
- Bug fix: Address non-deterministic behavior caused by inconsistent ordering of reads with identical start positions but different query names in sorted BAM files
- Bug fix: Fixed missing font in double arrow rendering

## [3.0.0]

- **Breaking change**: Modified how TRGT detects repeat motifs
  - Changed repeat segmentation algorithm to only match perfect STR motifs (imperfections in VNTR motifs are still allowed)
  - Implemented more flexible motif matching that allows better detection and visualization of repeat interruptions
- Made a number of improvements to TRGT plots
  - Replaced the size scale with the overall allele length
  - Implemented more accurate motif matching in waterfall plots
  - Made it possible to compress the plots vertically (useful for high-depth targeted data)
  - Made it easier to distinguish between mismatches and unsegmented regions
- Bug fix: Explicitly set reference from genome path for CRAM input, avoiding issues with moved/renamed references

## [2.1.0]

- **Important changes**:
  - Extended the use of optimized aligners to the plotting (TRVZ) functionality
  - Further optimized alignments, significantly improving overall performance for targeted datasets beyond the improvements introduced in 2.0, and cumulatively achieving ~2x speed-up for whole-genome sequencing data

## [2.0.0]

- **Important changes**:
  - Introduced fine-grained parallelization using work-stealing, significantly increasing CPU utilization, resulting in speed-ups of up to 200x for targeted datasets
  - All previous aligners used for genotyping have been replaced with significantly faster and more memory-efficient alternatives (WFA2-lib)
  - When using the targeted data preset, reads with higher repeat purity between flanks are now given priority
- TRGT now includes a default monospace font, ensuring consistent text rendering during plotting even on systems without fonts in standard locations.
- Added CLI option `--vcf-list` to TRGT merge, allowing users to specify a file containing a list of VCF files. This option is mutually exclusive with the existing `--vcf` flag.
- Introduced CLI option `--font-family` for specifying custom font families when plotting
- Bug fix: the VCF QUAL field is now always set with `.` rather than `0`
- Bug fix: fixed handling of AL tags when loading spanning BAM files modified by samtools
- Improved error handling and messaging during parsing of plotting inputs

## [1.5.1]

- Fixed an issue that prevented extraction of CpG methylation from BAM records containing multiple base modifications

## [1.5.0]

- Read clustering genotyper (`--genotyper cluster`) is now significantly faster at genotyping high-coverage repeat expansions; this may result in minor changes to consensus sequence and read assignment for highly mosaic repeats

## [1.4.1]

- Bug fix: corrected the type of the rq tag in BAM output

## [1.4.0]

- Parameters appropriate for targeted sequencing can now be set with `--preset targeted` option
- Waterfall plots no longer panic when there are no reads in a locus
- Algorithmic changes to `--genotyper cluster` allow fewer reads to be assigned to an allele; this may result in minor changes to consensus sequence and read assignment

## [1.3.0]

- Plotting code has been refactored as we prepare to revamp repeat visualizations
- The maximum number of reads per allele to plot can now be specified by `--max-allele-reads`
- Bug fix: repeat identifiers are now permitted to contain commas

## [1.2.0]

- `trgt merge`:
  - Multi-sample VCF Merging: Added support for merging TRGT VCFs with any number of samples, allowing updates to large, population-scale datasets with new samples
  - Synced contig indexing: Introduced support for VCFs with inconsistent contig orderings. Additionally the new `--contigs` flag allows specifying a comma-separated list of contigs to be merged
  - The reference genome is no longer required when merging TRGT VCFs from version 1.0.0 or later
  - Merging now skips and logs problematic loci by default. Use the `--quit-on-errors` flag to terminate on errors. Statistics are logged post-merge, including counts of failed and skipped TRs
- `trgt validate`
  - Always outputs statistics directly to stdout and stderr instead of logging them
- Bug fix: resolved issue with handling bgzip-compressed BED files

## [1.1.2]

- Bug fix: Prevent genotyping without reads
- Added the `--disable-bam-output` flag to `trgt genotype`, allowing users to disable BAMlet generation. **However, please note that BAMlets are still required for downstream tasks like trgt plot**

## [1.1.1]

- Bug fix: Read filtering logic no longer removes reads without RQ tags

## [1.1.0]

- Added a new subcommand `trgt merge`. This command merges VCF files generated by `trgt genotype` into a joint VCF file. **Works with VCFs generated by all versions of TRGT** (the resulting joint VCF will always be in the TRGT ≥v1.0.0 format which includes padding bases)
- Added subsampling of regions with ultra-high coverage (`>MAX_DEPTH * 3`, by default 750); implemented via reservoir sampling
- Fixed a cluster genotyper bug that occurred when only a single read covered a locus
- Added new logic for filtering non-HiFi reads: remove up to 3% of lower quality reads that do not match the expected repeat sequence

## [1.0.0]

- **Breaking change**: TRGT and TRVZ are now merged into a single binary. Users need to run subcommands `trgt genotype` and `trgt plot` for genotyping and visualization, respectively
- **Breaking change**:  A padding base is now automatically added to all genotyped allele sequences in the VCF file, ensuring better compliance with VCF standards and handling of zero-length alleles
- Added a new subcommand `trgt validate`. This command allows for validation of a repeat catalog against a given reference genome and reports statistics for any malformed entries
- Lower memory footprint: Better memory management significantly reduces memory usage with large repeat catalogs
- Updated error handling: Malformed entries are now logged as errors without terminating the program
- Added shorthand CLI options to simplify command usage

## [0.9.0]

- Add support for polyalanine repeats (by allowing characters `N` in the motif sequence)
- Fix a bug causing TRVZ to error out on polyalanine repeats

## [0.8.0]

- **Breaking change**: Motif spans and counts (`MS` and `MC` fields) and purity assessment (`AP` field) are now performed with an HMM-based algorithm for all repeats; expect some differences in results relative to the previous versions
- Allele purity of zero-length alleles are now reported as missing values in the VCFs
- The spanning.bam output file now carries over the QUAL values and mapping strand from the input reads
- Added an advanced flag `--output-flank-len` that controls the number of flanking bases reported in the spanning.bam files and shown in trvz plots
- A crash that may occur on BAMs where methylation was called twice has been fixed
- Optimizations to the `--genotyper=cluster` mode, including haploid genotyping of the X chromosome when `--karyotype` is set to `XY`

## [0.7.0]

- Read phasing information can now be used during repeat genotyping (via `HP` tags)
- Users can now define complex repeats by specifying motif sequences in the MOTIFS field and setting STRUC to <`locus_name`>
- The original MAPQ values in the input reads are now reported in the BAM output
- BAMlet sample name can now be provided using the `--sample-name` flag; if it not provided, it is extracted from the input BAM or file stem (addressing issue #18)

## [0.6.0]

- Add alignment CIGARs to spanning.bam reads
- Increase read extraction region
- Cluster genotyper reports confidence intervals
- Improved error handling of invalid input files (genome, catalog
    and reads)

## [0.5.0]

- The genotyper now uses information about SNPs adjacent to repeats
- BAM files now contain read-to-allele assignments
- Added support for gzip compressed repeat files
- Improved error handling and error messages

## [0.4.0]

- Added TRVZ tutorial
- Added sample karyotype parameter (`XX` or `XY`)
- Renamed VCF genotype field `ALCI` to `ALLR`
- Made genotyping algorithm changes to improve accuracy

## [0.3.4]

- Improved label spacing in TRVZ plots
