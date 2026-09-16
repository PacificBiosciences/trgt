# TRGT docker image

Builds the `trgt` binary from this fork (`bw2/trgt`) into a slim Debian image published as
`weisburd/trgt:latest` on Docker Hub. The digest of the last published image is recorded in
`sha256_dockerhub.txt`; pin that digest rather than `:latest` when a pipeline needs a fixed version.

## How it is published

`.github/workflows/build_docker_image.yml` runs after the `unit tests` workflow succeeds on a push to
`main`, builds the image from that tested commit, pushes it, and commits the updated
`sha256_dockerhub.txt` back to the repo. It can also be run by hand from the Actions tab
(https://github.com/bw2/trgt/actions/workflows/build_docker_image.yml), which skips the test gate.

It needs two repository secrets, set at https://github.com/bw2/trgt/settings/secrets/actions:

| Secret | Value |
| --- | --- |
| `DOCKERHUB_USERNAME` | Docker Hub account name |
| `DOCKERHUB_TOKEN` | Docker Hub access token with read/write scope |

## Building by hand

```bash
cd docker
make build          # docker build -f Dockerfile -t weisburd/trgt:latest ..
make push           # pushes and rewrites sha256_dockerhub.txt
```

The build context is the repository root, not this directory, because the image compiles the
checked-out source. `../.dockerignore` keeps the untracked `target/` and `benchmarks/` directories out
of that context, and deliberately keeps `.git` in it: `build.rs` runs `git describe --always`, and the
short commit sha it returns is what makes `trgt --version` report `5.1.0-<sha>` instead of a bare
`5.1.0`.

## Using it

There is no `ENTRYPOINT`, so the binary is invoked explicitly:

```bash
docker run --rm weisburd/trgt:latest trgt --version
docker run --rm -v "$PWD:/data" weisburd/trgt:latest \
    trgt genotype --genome /data/ref.fa --reads /data/sample.bam --repeats /data/catalog.bed --output-prefix /data/out
```
