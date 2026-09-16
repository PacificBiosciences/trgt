# Builds the TRGT binary from this repository (the bw2/trgt fork) into a slim runtime image.
#
# The build context is the repository ROOT, not this directory, because the source has to be copied
# in. Build it with `make build` from here, or directly:
#
#   docker build -f docker/Dockerfile -t weisburd/trgt:latest .
#
# See ../.dockerignore for what is kept out of that context (and why .git is kept in).

FROM rust:1.87-bookworm AS builder

# Native build dependencies of the two -sys crates:
#   hts-sys   compiles htslib with the curl/s3/gcs features enabled in Cargo.toml, so it needs zlib,
#             OpenSSL and libcurl headers. bzip2/lzma are NOT needed: rust-htslib is pulled with
#             default-features = false, and Cargo.lock accordingly has no bzip2-sys or lzma-sys.
#   wfa2-sys  builds the bundled WFA2-lib through cmake.
# Both run bindgen, which needs libclang at build time.
RUN apt-get update && apt-get install --no-install-recommends --yes \
    build-essential \
    cmake \
    pkg-config \
    libclang-dev \
    zlib1g-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .

# --locked builds exactly the dependency versions in the committed Cargo.lock and fails if the lock
# file is stale, rather than silently resolving newer ones.
RUN cargo build --release --locked \
    && strip target/release/trgt


FROM debian:bookworm-slim

# Runtime shared libraries for the htslib features compiled in above (reading remote BAM/CRAM over
# https, s3 and gcs), plus the CA bundle those TLS connections need.
RUN apt-get update && apt-get install --no-install-recommends --yes \
    ca-certificates \
    libcurl4 \
    libssl3 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/trgt /usr/local/bin/trgt

# Fail the build here rather than at run time if the runtime image is missing a shared library the
# binary was linked against. Also records the built version in the build log.
RUN trgt --version

# Commit this image was built from, passed by the workflow. Declared after the layers above so that
# changing it does not invalidate them.
ARG TRGT_COMMIT=unknown
LABEL org.opencontainers.image.revision="${TRGT_COMMIT}"
LABEL org.opencontainers.image.source="https://github.com/bw2/trgt"
LABEL org.opencontainers.image.description="TRGT built from the bw2/trgt fork"

# No ENTRYPOINT on purpose: these images get run both as `docker run <image> trgt genotype ...` and
# from workflow engines that invoke a shell inside the container, and an ENTRYPOINT gets in the way
# of the second. This matches how the str-analysis images behave.
CMD ["trgt", "--help"]
