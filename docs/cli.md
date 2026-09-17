# Command-line interface of TRGT

Commands:

- `genotype`
- `plot`
- `merge`
- `validate`
- `deepdive`

Options:

- `-v, --verbose` Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)
- `--color <COLOR>` Enable or disable color output in logging (always, auto, never); supports NO_COLOR/FORCE_COLOR environment variables [default: auto]
- `-h, --help` Print help
- `-V, --version` Print version

Most file input arguments accept remote URIs (`http://`, `https://`, `s3://`,
`gs://`). See [remote file support](remote_files.md) for details and
authentication tips.

## Environmental variables

- `TRGT_ENABLE_HTSLIB_LOGGING` Set this environment variable to re-enable HTSlib error messages, which are disabled by default in TRGT. This can be helpful for debugging unknown errors that might be related to HTSlib.
- `NO_COLOR` Set this environment variable to disable color output.
- `FORCE_COLOR` Set this environment variable to force color output.

For cloud storage and remote file access, additional environment variables are
available. See [remote file support](remote_files.md) for more details.

## Genotype command-line

Options:

- `-g, --genome <GENOME>` Path to the FASTA file containing reference genome. This must be the same reference genome as the one used for read alignment
- `-r, --reads <READS>` BAM file with aligned HiFi reads
- `--fail-reads <FAIL_READS>` Optional BAM file with fail reads; disables the default `rq` filter unless `--min-read-quality` is set
- `-b, --repeats <REPEATS>` BED file with reference coordinates and structure of tandem repeats
- `-o, --output-prefix <OUTPUT_PREFIX>` Prefix for output files. TRGT generates an unsorted VCF file (`<OUTPUT_PREFIX>.vcf.gz`) and unsorted BAM file with pieces of HiFi reads overlapping the repeats (`<OUTPUT_PREFIX>.spanning.bam`)
- `-k, --karyotype <KARYOTYPE>` Sample karyotype (XX, XY, or a file path) [default: XX]
- `-t, --threads <THREADS>` Number of threads [default: 1]
- `--preset <PRESET>` Parameter preset (wgs or targeted) [default: wgs]

Advanced:

- `--sample-name <SAMPLE_NAME>` Sample name. If not provided, the sample name is extracted from the input BAM or file stem
- `--genotyper <GENOTYPER>` Genotyping algorithm (size or cluster) [default: size; targeted preset default: cluster]
- `--flank-len <FLANK_LEN>` Minimum length of the flanking sequence [default: 250]
- `--output-flank-len <FLANK_LEN>` Length of flanking sequence to report on output [default: 50]
- `--disable-bam-output` Disable BAM output
- `--resume` Resume a previous run (see [resuming a run](#resuming-a-run)). Requires `--disable-bam-output`
- `--max-depth <MAX_DEPTH>` Maximum locus depth [default: 250; targeted preset default: 10000]
- `--fetcher-threads <THREADS>` Number of threads for querying input BAM files [default: half the number of analysis threads, max 8]

A locus with more than three times `--max-depth` reads keeps a sample of them,
and this fork of TRGT decides which ones from the locus itself: its contig, its
coordinates and its ID. Upstream TRGT decides it once for each group of nearby
loci it reads together, so the two produce different genotypes at loci deep
enough to be sampled. Deciding it per locus is what lets `--resume` (see
[resuming a run](#resuming-a-run)) finish a run with the genotypes an
uninterrupted run would have produced, and it also means the sample no longer
follows `--max-group-span` or `--max-group-size`.

Preset-dependent defaults:

- `--preset targeted` changes several defaults:
  - `--genotyper`: default becomes `cluster` (was `size`)
  - `--flank-len`: default becomes `200` (was `250`)
  - `--max-depth`: default becomes `10000` (was `250`)
  - Advanced (hidden):
    - `--aln-scoring`: default becomes `1,0,1` (was `2,5,1`)
    - `--min-flank-id-frac`: default becomes `0.8` (was `0.7`)
    - `--min-read-quality`: default becomes `-1.0` (was `0.98`)
- `--fail-reads` also changes `--min-read-quality` to `-1.0` unless `--min-read-quality` is specified explicitly.

### Resuming a run

A run that was interrupted, for example by a job time limit or a preempted
machine, can be continued with `--resume`. TRGT reads the output file the
interrupted run left behind, keeps the genotypes in it, and genotypes only the
loci of the catalog that are missing from it. Rerun the original command with
`--resume` added:

```bash
trgt genotype --genome reference.fasta \
              --repeats repeats.bed \
              --reads sample.bam \
              --output-prefix sample \
              --disable-bam-output \
              --resume
```

Points worth knowing:

- `--resume` requires `--disable-bam-output`, because TRGT cannot continue the
  `<OUTPUT_PREFIX>.spanning.bam` output. Resuming a run that wrote spanning
  reads would leave a BAM holding the reads of only the loci genotyped by the
  last run. A `<OUTPUT_PREFIX>.spanning.bam` that the interrupted run already
  wrote stays as that run left it, so it covers fewer loci than the resumed VCF
  beside it and `trgt plot` has no reads to draw for the rest. TRGT warns when it
  finds one. Delete it, or rerun without `--resume`, when the spanning reads are
  needed.
- Loci are matched by the `ID` field of the repeat catalog, so resume a run with
  the catalog it started with, and keep repeat IDs unique, as TRGT requires
  anyway (see [repeat files](repeat_files.md)). Before it touches the output
  file, a resumed run reads the catalog and checks it against the records it
  kept: each kept record has to match a catalog locus of the same ID covering the
  same region. It stops with an error naming the locus when the catalog has moved
  it, resized it, dropped it, or gives its ID to more loci than there are records
  to account for, since a record names its locus by ID alone and TRGT writes
  records in the order genotyping finishes them, so nothing would say which locus
  such a record belongs to. Loci the previous run never reached are not checked
  this way; a catalog changed only in those is resumed without complaint.
- TRGT refuses to resume an output file written by a different version of
  itself, since genotypes of two versions do not belong in one file. The
  `##trgtCommand` line of each run that contributed records is kept in the header
  of the resumed file.
- TRGT keeps the records of the previous run up to the first one it cannot read,
  and always genotypes the last of them again. A run killed while writing leaves
  a record that stops partway through, and a record whose text stops inside its
  last column still reads as a complete record, with the fields it no longer
  reaches read as missing. Genotyping that locus again costs one locus of work
  and keeps a half-written genotype out of the output. Records are compressed in
  blocks of up to 64 kB, so for a compressed output file the loci of the whole
  last block are genotyped again.
- TRGT refuses to resume an output file that holds a different sample, that was
  written against different reference contigs (compared against the contigs of
  the reads, since that is what the output file records), or that was genotyped
  with different settings or from different inputs. The `##trgtParams` header
  line of a run records its settings, along with a value standing for its reads
  and one standing for its reference, and the whole line has to match, so a
  resumed run cannot quietly finish a file with, say, a different `--karyotype`
  or a different `--reads` than the one it was started with. Neither value
  involves a path, so finding the same inputs somewhere else, as a run picked up
  on another machine does, is not a mismatch. Both stop short of the data itself,
  which would have to be read through to cover: the reads are stood for by the
  headers of their files, so two files holding different reads under the same
  header are not told apart, and the reference is stood for by the contigs of its
  index with their lengths, so a genome and a masked copy of it are not told
  apart either. Each refusal reports what differs and stops without touching the
  file. An output file that is unreadable altogether is reported and then
  overwritten, and the whole catalog is genotyped.
- A locus with more reads than `--max-depth` allows keeps a sample of them, and
  which reads it keeps is decided from the locus itself: where it lies and what
  it is called. A resumed run reads only the loci it still has to genotype, in
  groups that hold different neighbours than the groups of the run it continues,
  so a sample decided per group would have made those loci come out differently
  than an uninterrupted run would. Deciding it per locus is what keeps the two
  the same. It also means the reads a locus keeps no longer follow
  `--max-group-span` or `--max-group-size`. Note that this differs from upstream
  TRGT, which decides the sample once per group, so the two disagree on loci deep
  enough to be sampled, which is those above three times `--max-depth`.
- A run killed in the first moments of resuming, before it has written enough to
  fill an output block, loses the records it had just copied and genotypes those
  loci again. The loss is at most one block of records (a few hundred), and it is
  the whole checkpoint when the checkpoint is smaller than a block. Nothing is
  corrupted by it: the next resume picks up from what survived, and the output
  ends up the same as a run that was never interrupted.
- Resuming rewrites the records that are kept, so the output prefix needs free
  disk space for a second copy of the output file while the run starts up, and it
  reads the repeat catalog a second time to check it against those records before
  it replaces the file. For a catalog held remotely, that second read is a second
  download.

## Plot command-line

Options:

- `-g, --genome <GENOME>` Path to the FASTA file containing reference genome
- `-b, --repeats <REPEATS>` BED file with repeat coordinates and structure
- `-f, --vcf <VCF>` VCF file generated by TRGT
- `-r, --spanning-reads <SPANNING_READS>` BAM file with spanning reads generated by TRGT
- `-i, --repeat-id <REPEAT_ID>` ID of the repeat to visualize
- `-o, --image <OUTPUT_PATH>` Output image path that ends with extension .pdf, .svg, or .png

Plotting:

- `--plot-type <PLOT_TYPE>` Two types of plots can be generated: allele plots and waterfall plots. Allele plots show alignments of reads to each repeat allele. Waterfall plots display unaligned repeat sequences without aligning them to the (consensus) allele. Waterfall plots are especially useful for QC of repeat calls and for visualization of mosaic expansions [default: allele]
- `--squished` Horizontally compress the plot; useful for visualization of high-coverage targeted data
- `--show <SHOW>` Either motifs (`motifs`) or methylation (`meth`) is visualized [default: motifs]
- `--font-family <FONT>` Font family to use for text elements [default: Roboto Mono]

Advanced:

- `--flank-len <FLANK_LEN>` Length of flanking regions [default: 50]
- `--max-allele-reads <MAX_READS>` Max number of reads per allele to plot

## Merge command-line

The `merge` subcommand currently requires local files.

Options:

- `--vcf <VCF>` VCF files to merge
- `--vcf-list <VCF_LIST>` File containing paths of VCF files to merge (one per line)
- `-g, --genome <FASTA>` Path to the FASTA file containing reference genome
- `-o, --output <FILE>` Write output to a file [default: standard output]

Advanced:

- `-O, --output-type <OUTPUT_TYPE>` Output type: u|b|v|z, u/b: un/compressed BCF, v/z: un/compressed VCF
- `--skip-n <SKIP_N>` Skip the first N records
- `--process-n <PROCESS_N>` Only process N records
- `--print-header` Print only the merged header and exit
- `--force-single` Run even if there is only one file on input
- `--no-version` Do not append version and command line to the header
- `--quit-on-errors` Quit immediately on errors during merging
- `--contig <CONTIG>` Process only the specified contigs (comma-separated list)
- `-t, --threads <THREADS>` Number of threads for (de)compressing input/output VCF files (shared across readers and writer) [default: 2]
- `--no-index` Stream VCFs without loading their indexes (contig order must match across inputs)
- `-W, --write-index` Write index for the output compressed VCF/BCF file

## Validate command-line

Options:

- `-g, --genome <FASTA>` Path to the FASTA file containing reference genome
- `-b, --repeats <REPEATS>` BED file with repeat coordinates and structure

Advanced:

- `--flank-len <FLANK_LEN>` Length of flanking regions [default: 50]

## Deepdive command-line

Options:

- `-g, --genome <FASTA>` Path to reference genome FASTA
- `-f, --vcf <VCF>` VCF file generated by trgt genotype
- `-b, --repeats <REPEATS>` BED file with repeat coordinates
- `-r, --spanning-reads <SPANNING_READS>` BAM file with spanning reads generated by trgt genotype
- `-i, --repeat-id <REPEAT_ID>` ID of the repeat to realign
- `-o, --output-prefix <OUTPUT_PREFIX>` Prefix for output files (.fasta, .bed and .bam)
