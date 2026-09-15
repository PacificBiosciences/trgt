# Remote file support

TRGT can read inputs directly from remotely hosted files. Anywhere the CLI
accepts a local path you can usually provide a URI instead. Outputs are always
written to the local filesystem.

## Files that can be remote

- Reference genomes (`--genome`)
- Read alignments (`--reads`, `--fail-reads`, `--spanning-reads`)
- Repeat catalogs (`--repeats`)
- VCF files for plotting and deepdive (`--vcf`)

The `merge` subcommand currently requires local files.

## Supported URI schemes

TRGT accepts `http://`, `https://`, `s3://`, and `gs://` locations in addition
to normal filesystem paths. At runtime TRGT performs validation and
logs warnings for common misconfigurations (for example missing credentials).

## Performance and reliability

Remote file support is still experimental. Performance may be slower when working with large
FASTA, BAM/CRAM, VCF, and/or BED files, especially over long-distance networks.

## HTTPS certificates and CA bundles

TRGT uses `libcurl` to access remote files. When connecting to HTTPS endpoints
or cloud gateways, `libcurl` must verify the server's security certificate using
a Certificate Authority (CA) bundle. TRGT automatically tries to find a system
CA bundle and configure `libcurl` to use it.

- To use a specific CA file, set `CURL_CA_BUNDLE=/path/to/cacert.pem`.
- To disable automatic detection, set `TRGT_DISABLE_CA_AUTODETECT=1`.

If no bundle can be found, you might see an error like: `problem with the SSL CA cert (60/77)`. 
To fix this, install a system CA bundle or download `cacert.pem` from the
[cURL project](https://curl.se/docs/caextract.html), then point `CURL_CA_BUNDLE`
to that file.

## Cloud authentication

When accessing cloud object stores, TRGT relies on the same authentication
methods as aws, gsutil, and other standard tools.

- **Amazon S3** (`s3://`): Use AWS settings such as `AWS_PROFILE`,
  `AWS_ACCESS_KEY_ID`/`AWS_SECRET_ACCESS_KEY`, or a shared credential file. For
  public buckets, set `AWS_SHARED_CREDENTIALS_FILE=/dev/null` to disable credential lookup and enable anonymous access.
- **Google Cloud Storage** (`gs://`): Set `GOOGLE_APPLICATION_CREDENTIALS` to a
  service account JSON key file or authenticate interactively using `gcloud auth application-default login`. 
  If neither is configured, TRGT will warn on runtime.

Access still depends on the bucket’s permissions and whether the system running
TRGT can reach the endpoint.

### Environment variables

| Provider | Purpose | Variable(s) |
|----------|---------|-------------|
| HTTPS | Override CA trust store | `CURL_CA_BUNDLE` |
| HTTPS | Disable CA auto-detection | `TRGT_DISABLE_CA_AUTODETECT` |
| Amazon S3 | Access keys / profiles | `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_PROFILE` |
| Amazon S3 | Anonymous access to public buckets | `AWS_SHARED_CREDENTIALS_FILE=/dev/null` |
| Google Cloud Storage | Service account credentials | `GOOGLE_APPLICATION_CREDENTIALS` |

## Quick start example

The command below runs `trgt genotype` while streaming inputs from three remote
providers: an S3-hosted reference genome, a Google Cloud Storage BAM file, and
an HTTPS repeat catalog. The S3 bucket is public, so the command sets
`AWS_SHARED_CREDENTIALS_FILE=/dev/null` to disable local credentials and enable
anonymous access. If your environment requires authentication, remove this
override or point it to your credentials file.

```bash
AWS_SHARED_CREDENTIALS_FILE=/dev/null \
trgt genotype \
  --genome 's3://pacbio-hifi-human-wgs-reference/dataset/GRCh38/human_GRCh38_no_alt_analysis_set.fasta' \
  --reads 'gs://deepvariant/pacbio-case-study-testdata/HG003.pfda_challenge.35x.grch38.bam' \
  --repeats 'https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.hg38.TRGT.bed' \
  --output-prefix tmp \
  --preset targeted \
  --threads 4 \
  -v
```

## Troubleshooting

- Remote file access fails: set `TRGT_ENABLE_HTSLIB_LOGGING=1` to enable detailed error messages from the underlying HTSlib library. This can help identify permission errors, missing files, or authentication problems.
- `SSL peer certificate or SSH remote key was not OK` or `problem with the SSL CA cert`: provide a valid `CURL_CA_BUNDLE` or disable auto-detection if using a custom trust store.
- `Failed to create BAM reader from URI`: Check that the URI is correct, the object exists, and your credentials allow read access.
- `The 'GOOGLE_APPLICATION_CREDENTIALS' environment variable is not set`: Authenticate using `gcloud auth application-default login` or point to a service account key file.
- Slow downloads or long runtime: High-latency networks can slow-down streaming. Copy large inputs locally for repeated analyses.

If you encounter issues with remote access, please report them so this feature
can continue improving.
