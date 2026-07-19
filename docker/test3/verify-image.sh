#!/usr/bin/env bash
set -Eeuo pipefail

IMAGE_NAME="${IMAGE_NAME:-spatialsnake-test3}"
IMAGE_TAG="${IMAGE_TAG:-0.0.2-test3-20260719}"
IMAGE_REF="${IMAGE_NAME}:${IMAGE_TAG}"

docker image inspect "${IMAGE_REF}" >/dev/null

docker run --rm \
    --entrypoint python \
    "${IMAGE_REF}" \
    /usr/local/lib/spatialsnake/smoke_test.py

docker run --rm "${IMAGE_REF}" --version

MOUNT_TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/spatialsnake-bind-test.XXXXXX")"
trap 'rm -rf "${MOUNT_TEST_DIR}"' EXIT

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --mount "type=bind,src=${MOUNT_TEST_DIR},dst=/work" \
    --entrypoint bash \
    "${IMAGE_REF}" \
    -lc 'printf "container bind mount OK\n" > /work/result.txt'

grep -qx 'container bind mount OK' "${MOUNT_TEST_DIR}/result.txt"

if [[ "${VERIFY_GPU:-0}" == "1" ]]; then
    docker run --rm \
        --gpus all \
        --entrypoint python \
        "${IMAGE_REF}" \
        -c 'import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name(0))'
fi

echo "Verified ${IMAGE_REF}"
