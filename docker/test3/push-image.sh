#!/usr/bin/env bash
set -Eeuo pipefail

SOURCE_IMAGE_NAME="${IMAGE_NAME:-spatialsnake-test3}"
IMAGE_TAG="${IMAGE_TAG:-0.0.2-test3-20260719}"
SOURCE_IMAGE="${SOURCE_IMAGE_NAME}:${IMAGE_TAG}"

REGISTRY="${REGISTRY:-docker.io}"
REGISTRY_NAMESPACE="${REGISTRY_NAMESPACE:-}"
REGISTRY_REPOSITORY="${REGISTRY_REPOSITORY:-spatialsnake-test3}"

if [[ -z "${REGISTRY_NAMESPACE}" ]]; then
    echo "Error: set REGISTRY_NAMESPACE to your Docker Hub user or registry namespace." >&2
    exit 1
fi

if [[ "${REGISTRY}" == "docker.io" ]]; then
    TARGET_IMAGE="${REGISTRY_NAMESPACE}/${REGISTRY_REPOSITORY}:${IMAGE_TAG}"
else
    TARGET_IMAGE="${REGISTRY}/${REGISTRY_NAMESPACE}/${REGISTRY_REPOSITORY}:${IMAGE_TAG}"
fi

docker image inspect "${SOURCE_IMAGE}" >/dev/null
docker image tag "${SOURCE_IMAGE}" "${TARGET_IMAGE}"
docker image push "${TARGET_IMAGE}"

echo "Pushed ${TARGET_IMAGE}"
