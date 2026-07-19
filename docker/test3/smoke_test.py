#!/usr/bin/env python
from __future__ import annotations

import importlib
from pathlib import Path
import subprocess
import sys


PYTHON_MODULES = [
    "spatialsnake",
    "snakemake",
    "spatialdata",
    "spatialdata_io",
    "spatialdata_plot",
    "scanpy",
    "squidpy",
    "cell2location",
    "scvi",
    "torch",
    "cellphonedb",
    "cellcharter",
    "banksy",
    "liana",
    "pyscenic",
    "pydeseq2",
    "bbknn",
    "geosketch",
]

R_PACKAGES = [
    "optparse",
    "Seurat",
    "spacexr",
    "Matrix",
    "dplyr",
    "ggplot2",
    "tibble",
    "pheatmap",
    "schard",
    "CellChat",
    "jsonlite",
    "future",
    "patchwork",
    "ComplexHeatmap",
    "AnnotationDbi",
    "clusterProfiler",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "ggsankey",
    "cowplot",
]

ENV_BIN = Path(sys.executable).resolve().parent


def check_python() -> None:
    for module_name in PYTHON_MODULES:
        module = importlib.import_module(module_name)
        version = getattr(module, "__version__", "unknown")
        print(f"PYTHON OK\t{module_name}\t{version}")


def check_r() -> None:
    package_vector = ", ".join(f'"{package}"' for package in R_PACKAGES)
    expression = (
        f"pkgs <- c({package_vector}); "
        "for (p in pkgs) { "
        "if (!requireNamespace(p, quietly=TRUE)) stop(paste('missing R package:', p)); "
        "cat('R OK\\t', p, '\\t', as.character(packageVersion(p)), '\\n', sep='') "
        "}"
    )
    subprocess.run([str(ENV_BIN / "Rscript"), "-e", expression], check=True)


def check_cli() -> None:
    subprocess.run([str(ENV_BIN / "spatialsnake"), "--version"], check=True)
    subprocess.run([str(ENV_BIN / "snakemake"), "--version"], check=True)


def main() -> int:
    print(sys.version)
    check_python()
    check_r()
    check_cli()
    print("Spatialsnake test3 runtime smoke test passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
