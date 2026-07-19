#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCK_DIR="${SCRIPT_DIR}/locks"
ARTIFACT_DIR="${SCRIPT_DIR}/artifacts"
ENV_NAME="${CONDA_ENV_NAME:-test3}"
ENV_PREFIX="${CONDA_ENV_PREFIX:-}"
ARCHIVE="${ARTIFACT_DIR}/test3-linux-64.tar.gz"
COMPRESS_LEVEL="${COMPRESS_LEVEL:-1}"

CONDA_BIN="${CONDA_EXE:-$(command -v conda || true)}"
CONDA_PACK_BIN="${CONDA_PACK_EXE:-$(command -v conda-pack || true)}"
CONDA_PACK_PYTHON="${CONDA_PACK_PYTHON:-}"

if [[ -z "${CONDA_BIN}" ]]; then
    echo "Error: conda was not found in PATH." >&2
    exit 1
fi

if [[ -z "${ENV_PREFIX}" ]]; then
    ENV_PREFIX="$("${CONDA_BIN}" env list | awk -v env_name="${ENV_NAME}" '$1 == env_name {print $NF; exit}')"
fi

if [[ -z "${ENV_PREFIX}" || ! -d "${ENV_PREFIX}/conda-meta" ]]; then
    echo "Error: conda environment '${ENV_NAME}' was not found. Set CONDA_ENV_PREFIX explicitly." >&2
    exit 1
fi

if [[ -z "${CONDA_PACK_BIN}" && -z "${CONDA_PACK_PYTHON}" ]]; then
    cat >&2 <<'EOF'
Error: conda-pack was not found.
Install it outside the test3 environment, for example:
  conda install -n base -c conda-forge conda-pack
Or set CONDA_PACK_EXE / CONDA_PACK_PYTHON explicitly.
EOF
    exit 1
fi

if [[ -z "${CONDA_PACK_PYTHON}" ]]; then
    CONDA_PACK_PYTHON="$(sed -n '1s/^#!//p' "${CONDA_PACK_BIN}")"
fi

if [[ ! -x "${CONDA_PACK_PYTHON}" ]]; then
    echo "Error: unable to resolve the Python interpreter used by conda-pack." >&2
    echo "Set CONDA_PACK_PYTHON explicitly." >&2
    exit 1
fi

mkdir -p "${LOCK_DIR}" "${ARTIFACT_DIR}"

echo "Exporting conda metadata from ${ENV_PREFIX}"
"${CONDA_BIN}" export --prefix "${ENV_PREFIX}" \
    --format=environment-yaml \
    --file="${LOCK_DIR}/environment-linux-64.yml"
"${CONDA_BIN}" export --prefix "${ENV_PREFIX}" \
    --from-history \
    --format=environment-yaml \
    --file="${LOCK_DIR}/environment-from-history.yml"
if ! "${CONDA_BIN}" export --prefix "${ENV_PREFIX}" \
    --format=explicit \
    --file="${LOCK_DIR}/explicit-linux-64.txt"; then
    echo "Conda refused a mixed conda/pip explicit export; exporting the conda package subset." >&2
    "${CONDA_BIN}" list --prefix "${ENV_PREFIX}" --explicit \
        > "${LOCK_DIR}/explicit-linux-64.txt"
fi

"${ENV_PREFIX}/bin/python" -m pip freeze --all \
    > "${LOCK_DIR}/pip-freeze.txt"

"${ENV_PREFIX}/bin/Rscript" - "${LOCK_DIR}/r-packages.tsv" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
output <- args[[1]]
packages <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
fields <- c(
  "Package", "Version", "LibPath", "Repository",
  "RemoteType", "RemoteHost", "RemoteUsername", "RemoteRepo",
  "RemoteRef", "RemoteSha"
)
for (field in setdiff(fields, colnames(packages))) {
  packages[[field]] <- NA_character_
}
write.table(
  packages[, fields],
  file = output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
RS

{
    echo "exported_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "source_prefix=${ENV_PREFIX}"
    echo "source_os=$(uname -s)"
    echo "source_arch=$(uname -m)"
    echo "source_kernel=$(uname -r)"
    echo "source_glibc=$(getconf GNU_LIBC_VERSION 2>/dev/null || true)"
    echo "conda_version=$("${CONDA_BIN}" --version)"
    echo "python_version=$("${ENV_PREFIX}/bin/python" --version 2>&1)"
    echo "r_version=$("${ENV_PREFIX}/bin/Rscript" -e 'cat(R.version.string)' 2>&1)"
    echo "snakemake_version=$("${ENV_PREFIX}/bin/snakemake" --version 2>&1)"
} > "${LOCK_DIR}/system-info.txt"

echo "Packing ${ENV_PREFIX} into ${ARCHIVE}"
"${CONDA_PACK_PYTHON}" "${SCRIPT_DIR}/conda_pack_compat.py" \
    --prefix "${ENV_PREFIX}" \
    --output "${ARCHIVE}" \
    --compress-level "${COMPRESS_LEVEL}"

(
    cd "${SCRIPT_DIR}"
    {
        sha256sum artifacts/test3-linux-64.tar.gz
        find locks -maxdepth 1 -type f ! -name SHA256SUMS -print0 \
            | sort -z \
            | xargs -0 sha256sum
    } > locks/SHA256SUMS
)

echo "Created ${ARCHIVE}"
echo "Archive SHA-256: $(sha256sum "${ARCHIVE}" | awk '{print $1}')"
