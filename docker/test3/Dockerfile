# syntax=docker/dockerfile:1.7

ARG BASE_IMAGE=ubuntu:22.04@sha256:4f838adc7181d9039ac795a7d0aba05a9bd9ecd480d294483169c5def983b64d
FROM ${BASE_IMAGE}

ARG APP_UID=1000
ARG APP_GID=1000

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates tini \
    && rm -rf /var/lib/apt/lists/*

RUN groupadd --gid "${APP_GID}" spatialsnake \
    && useradd --uid "${APP_UID}" --gid "${APP_GID}" --create-home --no-log-init spatialsnake \
    && mkdir -p /opt/conda/envs/test3 /opt/spatialsnake-manifest /work \
    && chown -R spatialsnake:spatialsnake /work /home/spatialsnake

RUN --mount=type=bind,source=docker/test3/artifacts/test3-linux-64.tar.gz,target=/tmp/test3-linux-64.tar.gz,readonly \
    tar -xzf /tmp/test3-linux-64.tar.gz -C /opt/conda/envs/test3 \
    && PATH="/opt/conda/envs/test3/bin:${PATH}" /opt/conda/envs/test3/bin/conda-unpack

ENV PATH="/opt/conda/envs/test3/bin:${PATH}" \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/spatialsnake-matplotlib \
    NUMBA_CACHE_DIR=/tmp/spatialsnake-numba \
    XDG_CACHE_HOME=/tmp/spatialsnake-cache \
    HOME=/home/spatialsnake \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

WORKDIR /opt/spatialsnake-src

COPY setup.py MANIFEST.in README.md requirements.txt requirements-extended.txt ./
COPY spatialsnake ./spatialsnake
COPY docker/test3/locks /opt/spatialsnake-manifest
COPY docker/test3/smoke_test.py /usr/local/lib/spatialsnake/smoke_test.py

RUN python -m pip install --no-deps --no-build-isolation . \
    && python /usr/local/lib/spatialsnake/smoke_test.py \
    && rm -rf /tmp/spatialsnake-cache /tmp/spatialsnake-matplotlib /tmp/spatialsnake-numba

ARG IMAGE_VERSION=local
ARG BUILD_DATE=unknown
ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="Spatialsnake test3 runtime" \
      org.opencontainers.image.description="Spatialsnake with the packed test3 Python/R conda environment" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.created="${BUILD_DATE}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.source="https://github.com/l-zh007/spatialsnake" \
      org.opencontainers.image.licenses="MIT"

USER spatialsnake
WORKDIR /work

ENTRYPOINT ["/usr/bin/tini", "--", "spatialsnake"]
CMD ["--help"]
