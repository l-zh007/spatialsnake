#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ARCHIVE="${SCRIPT_DIR}/artifacts/test3-linux-64.tar.gz"

IMAGE_NAME="${IMAGE_NAME:-spatialsnake-test3}"
IMAGE_TAG="${IMAGE_TAG:-0.0.2-test3-20260719}"
IMAGE_REF="${IMAGE_NAME}:${IMAGE_TAG}"
APP_UID="${APP_UID:-1000}"
APP_GID="${APP_GID:-1000}"

if [[ ! -s "${ARCHIVE}" ]]; then
    echo "Error: ${ARCHIVE} is missing. Run docker/test3/export-test3.sh first." >&2
    exit 1
fi

if [[ -f "${SCRIPT_DIR}/locks/SHA256SUMS" ]]; then
    (
        cd "${SCRIPT_DIR}"
        sha256sum --check --ignore-missing locks/SHA256SUMS
    )
fi

BUILD_DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
VCS_REF="$(git -C "${PROJECT_ROOT}" rev-parse --short=12 HEAD 2>/dev/null || echo unknown)"

docker build \
    --pull \
    --file "${SCRIPT_DIR}/Dockerfile" \
    --tag "${IMAGE_REF}" \
    --build-arg "IMAGE_VERSION=${IMAGE_TAG}" \
    --build-arg "BUILD_DATE=${BUILD_DATE}" \
    --build-arg "VCS_REF=${VCS_REF}" \
    --build-arg "APP_UID=${APP_UID}" \
    --build-arg "APP_GID=${APP_GID}" \
    "${PROJECT_ROOT}"

echo "Built ${IMAGE_REF}"
