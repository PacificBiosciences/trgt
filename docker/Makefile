# Tag published by .github/workflows/build_docker_image.yml, and by `make push` from a laptop.
DOCKERHUB_TAG = weisburd/trgt:latest
DOCKERHUB_BASE = $(subst :latest,,$(DOCKERHUB_TAG))

# Commit baked into the image label. CI passes the commit it built; a plain `make build` picks up
# whatever is checked out here.
TRGT_COMMIT ?= $(shell git rev-parse HEAD 2>/dev/null || echo unknown)

all: build push

# The build context is the repository root (..), since the Dockerfile compiles this checkout.
build:
	docker build -f Dockerfile --build-arg TRGT_COMMIT=$(TRGT_COMMIT) -t $(DOCKERHUB_TAG) ..

push:
	docker push $(DOCKERHUB_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 | sed 's|^|$(DOCKERHUB_BASE)@|' > sha256_dockerhub.txt

.PHONY: all build push
