#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARCHIVE="${SCRIPT_DIR}/artifacts/test3-linux-64.tar.gz"
CHUNK_SIZE="${CHUNK_SIZE:-1900M}"
PREFIX="${ARCHIVE}.part-"

if [[ ! -s "${ARCHIVE}" ]]; then
    echo "Error: ${ARCHIVE} is missing. Run docker/test3/export-test3.sh first." >&2
    exit 1
fi

(
    cd "${SCRIPT_DIR}"
    sha256sum --check --ignore-missing locks/SHA256SUMS
)

rm -f "${PREFIX}"*
split -b "${CHUNK_SIZE}" -d -a 3 "${ARCHIVE}" "${PREFIX}"
sha256sum "${PREFIX}"* > "${ARCHIVE}.parts.sha256"

echo "Created chunks:"
ls -lh "${PREFIX}"*
echo "Chunk checksums: ${ARCHIVE}.parts.sha256"
