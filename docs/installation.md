# Installing Spatialsnake

These instructions target Linux and a fresh Conda environment. Activate the environment before using `python -m pip`. The PyPI route below uses Spatialsnake 0.0.6 with Python 3.12.11, R 4.4.0 and Snakemake 9.8.1.

## 1. PyPI: base installation

Create the environment, then install its R dependencies **before** installing Python packages with pip. Mamba may be substituted for Conda. Keep the Python/R/Snakemake constraints on the Bioconductor installation to avoid an unconstrained environment change.

```bash
conda create -n spatialsnake_env --override-channels -c conda-forge -c bioconda \
  python=3.12.11 snakemake-minimal=9.8.1 r-base=4.4.0 -y
conda activate spatialsnake_env

conda install --override-channels -c conda-forge -c bioconda \
  r-optparse r-tidyverse r-future r-jsonlite r-rcolorbrewer r-patchwork \
  r-cowplot r-pheatmap r-seurat r-remotes r-biocmanager r-presto r-nmf \
  r-circlize r-ggpubr -y

conda install --override-channels -c conda-forge -c bioconda \
  python=3.12.11 r-base=4.4.0 snakemake-minimal=9.8.1 \
  bioconductor-annotationdbi bioconductor-complexheatmap \
  bioconductor-clusterprofiler bioconductor-edger \
  bioconductor-org.hs.eg.db bioconductor-org.mm.eg.db \
  bioconductor-rhdf5 bioconductor-biocneighbors -y

python -m pip install spatialsnake==0.0.6
spatialsnake --version
python -m pip check
```

The R layer includes binary dependencies useful for later optional modules, avoiding repeated compilation. The base Python package and this R layer provide core analysis dependencies. Core analysis does not require PyTorch or a GPU. The repository's existing `environment.yml` contains only the initial Python/R/Snakemake environment; it does not replace the complete dependency steps above.

### Optional analyses: PyPI installations

Install **one** PyTorch build before installing the extended Python dependencies.

**CPU:**

```bash
python -m pip install --index-url https://download.pytorch.org/whl/cpu "torch==2.8.0"
python -m pip install "spatialsnake[extended]==0.0.6"
spatialsnake install-packages
python -m pip check
```

**Alternatively, an NVIDIA GPU compatible with the CUDA 12.8 PyTorch build:**

```bash
python -m pip install --index-url https://download.pytorch.org/whl/cu128 "torch==2.8.0"
python -m pip install "spatialsnake[extended]==0.0.6"
spatialsnake install-packages
python -m pip check
```

Do not run both backend selections. Installing the extended extra adds the Python modules; `install-packages` completes the additional BANKSY/R dependencies. In a PyPI/source installation, `install-packages --extended` does **not** replace `pip install "spatialsnake[extended]"`.

## 2. Bioconda: alternative package route

Use a fresh environment and the version available in the [Bioconda package listing](https://anaconda.org/bioconda/spatialsnake). This route installs the package through Conda rather than pip:

```bash
conda create -n spatialsnake_env --override-channels -c conda-forge -c bioconda \
  spatialsnake python=3.12.11 r-base=4.4.3 snakemake-minimal=9.8.1 -y
conda activate spatialsnake_env
spatialsnake install-packages
spatialsnake --version
python -m pip check
```

The documented Bioconda combination uses R 4.4.3, whereas the PyPI combination above uses R 4.4.0. These are installation-route constraints, not an analytical requirement to use different R versions. Record the version actually installed for reproducibility.

For optional modules, the **Bioconda** installer can select the PyTorch build and install the extended Python/R packages:

```bash
spatialsnake install-packages --extended --torch-backend=auto
python -m pip check
```

`auto` reuses an installed PyTorch build that satisfies the required version; it does not automatically download a build. If PyTorch is not installed, choose and install a CPU or CUDA build as shown in section 1, or use **one** of the explicit backend commands below:

```bash
spatialsnake install-packages --extended --torch-backend=cpu
```

```bash
spatialsnake install-packages --extended --torch-backend=cu128
```

`--extended` includes the Bioconda base-completion step. To inspect the planned installation without executing it, add `--dry-run`.

## 3. Source installation: development or a specified revision

Use the initial environment and R dependency steps from section 1, then install the checked-out source instead of the PyPI package:

```bash
git clone https://github.com/l-zh007/spatialsnake.git
cd spatialsnake
# For a reproducible checkout, select the release tag or commit for your analysis:
# git checkout <release-tag-or-commit>
python -m pip install .
spatialsnake --version
```

The default branch is a development checkout and can differ from a published package even when its version string has not changed. Record its commit and any local modifications. Do not assume a release tag exists without checking the repository.

For optional source dependencies, choose the CPU or CUDA PyTorch command from section 1, then run:

```bash
python -m pip install ".[extended]"
spatialsnake install-packages
python -m pip check
```

## 4. Which dependencies do I need?

| Analysis scope | Installation | GPU requirement |
|---|---|---|
| Core ingestion, preprocessing, clustering, marker/annotation support and common object utilities | Base package and the R dependency layer | Not required |
| cell2location and CellCharter | Extended Python packages and their dependencies | CPU or compatible GPU; GPU can accelerate model training |
| CellPhoneDB, LIANA, pySCENIC, Squidpy and optional viewer tools | Extended Python packages | No general GPU requirement; requirements depend on the selected method |
| BANKSY, RCTD, CellChat and Seurat export dependencies | Complete the corresponding installer steps above | No general GPU requirement |

Installation-time `--torch-backend` selects a PyTorch distribution. Analysis-time `--device=auto`, where supported, selects CUDA when available and CPU otherwise; it does not install a backend. Some R packages are obtained from GitHub and need network access. The installer detects installed dependencies and can be rerun after an interrupted download.

For installation problems, report the installation route, `spatialsnake --version`, Python/R versions, the command and its complete error log in a [GitHub issue](https://github.com/l-zh007/spatialsnake/issues).
