# Merging TRGT VCFs

The `trgt merge` command can be used to merge TRGT VCFs into a single multi-sample
VCF.

## Input requirements

By default, input TRGT VCFs must be sorted, BGZF-compressed, and indexed (`.tbi` or `.csi`).
All VCFs must be genotyped with the same repeat catalog and accessible on the
local filesystem. If any input VCF predates TRGT 1.0 the merge will require a
reference genome file to be added with `--genome` such that the missing padding
base can be reintroduced. Sample IDs must be unique across files, duplicates
cause the merge to abort and should be renamed before retrying.

## Output formatting

Without `--output` set, results will stream to standard output in uncompressed
VCF format. Set `--output` to write to a file. The file extension (`.vcf`,
`.vcf.gz`, `.bcf`, `.bcf.gz`) will automatically determine both format and
compression level unless you override it using `--output-type`. Use
`--print-header` to inspect the merged header without merging variants, this is
useful to verify sample order or merged metadata.

## How records are merged

Records are traversed contig by contig in reference order, and you can constrain
work to specific contigs with `--contig chr20,chrX`. For each site the reference
alleles must match across all inputs, otherwise the site fails with a clear
error. By default the merge logs the error and skips the record, add
`--quit-on-errors` to terminate immediately. FORMAT fields (`AL`, `ALLR`, `SD`,
`MC`, `MS`, `AP`, `AM`) are merged sample-by-sample, with missing samples filled
using the correct sentinel values so downstream tooling sees absence
explicitly. Headers are merged once up-front, the command rewrites historic
field definitions when needed (for example, collapsing `ALCI` into the current
`ALLR` representation) and injects version metadata unless `--no-version` is
set.

## Example usage

Merge three TRGT VCFs into one compressed VCF:

```bash
trgt merge \
  --vcf results/sample1.trgt.vcf.gz results/sample2.trgt.vcf.gz results/sample3.trgt.vcf.gz \
  --output cohort.trgt.vcf.gz
```

If the results folders contains these three VCF files you can also merge using auto-expansions,
while also writing an index file after merging and not loading the indices of the input VCFs:

```bash
trgt merge \
  --vcf results/*.vcf.gz \
  --output cohort.trgt.vcf.gz \ 
  --write-index \
  --no-index
```

Running the merge from a manifest file, while also limiting processing to chromosome 12:

```bash
trgt merge \
  --vcf-list cohort_vcfs.txt \
  --contig chr12 \
  --output cohort.chr12.trgt.vcf.gz
```

Here `cohort_vcfs.txt` is a plain text file with one VCF path per line. In this
file leading and trailing whitespace is ignored, blank lines are allowed, and
lines starting with `#` are treated as comments. Both relative and absolute
paths are accepted as long as they exist when the merge runs.

## Advanced usage

- `--threads N` sets the shared thread pool for decompressing inputs and
  compressing outputs (defaults to 2). Depending on the underlying disk device used
  this can significantly improve merging speed, especially when reading or writing
  compressed VCF/BCF.
- `--no-index` streams VCFs without `.tbi`/`.csi` files. Only use this when all
  inputs have identical contig order and are already position-sorted, otherwise
  the merge aborts with a sorting error. Streaming avoids loading index files
  into memory, which can cut RAM use significantly for large cohorts with many VCFs.

## Troubleshooting

- If an error is reported with "Reference alleles do not match", it indicates that
  the input VCFs used different TR catalogs for genotyping.
- High memory use: try `--no-index` to avoid loading `.tbi`/`.csi` into RAM, and
  consider hierarchical merging instead of a single monolithic run. For example: merge `a` + `b`
  -> `ab`, merge `c` + `d` -> `cd`, then merge `ab` + `cd` -> `abcd`.
