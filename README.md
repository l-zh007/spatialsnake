<p align="center"> 
  <img src="spatialsnake-logo.png" alt="spatialsnake logo" width="760">
</p>

<p align="center">
  <img
    width="820"
    src="https://readme-typing-svg.demolab.com?font=Fira+Code&weight=500&size=16&duration=3200&pause=1100&color=4C78A8&center=true&vCenter=true&width=820&lines=Unified+SpatialData-based+workflows+across+spatial+platforms;Decision-guided+modules+from+data+ingestion+to+annotation;Reproducible+single-sample+and+cross-sample+analysis;Modular+tools+for+integration%2C+transformation%2C+and+comparison"
    alt="Spatialsnake workflow features"
  />
</p>

<p align="center">
  <a href="https://spatialsnake-tutorial.readthedocs.io/en/latest/">
    <img src="https://img.shields.io/badge/docs-Read%20the%20Documentation-blue" alt="Documentation">
  </a>
  <a href="https://pypi.org/project/spatialsnake/">
    <img src="https://img.shields.io/pypi/v/spatialsnake" alt="PyPI">
  </a>
  <a href="https://anaconda.org/bioconda/spatialsnake">
    <img src="https://anaconda.org/bioconda/spatialsnake/badges/version.svg" alt="Bioconda version">
  </a>
  <a href="https://anaconda.org/bioconda/spatialsnake">
    <img src="https://anaconda.org/bioconda/spatialsnake/badges/platforms.svg" alt="Supported platforms">
  </a>
  <a href="https://anaconda.org/bioconda/spatialsnake">
    <img src="https://anaconda.org/bioconda/spatialsnake/badges/license.svg" alt="License">
  </a>
  <a href="https://anaconda.org/bioconda/spatialsnake">
    <img src="https://anaconda.org/bioconda/spatialsnake/badges/downloads.svg" alt="Bioconda downloads">
  </a>
</p>

<p align="center">
  A user-friendly command-line application for spatial transcriptomics powered by <code>SpatialData</code> and Snakemake.
</p>

`spatialsnake` is a command-line application for automated `spatial transcriptomics` analysis. Implemented in `Python` on top of the `scverse` ecosystem, it uses SpatialData to convert datasets from multiple spatial transcriptomics platforms into a unified `zarr`-based object format. This design supports a consistent workflow spanning data ingestion, preprocessing, clustering, annotation, and downstream analysis through a command-line interface with workflow-based parameter control.

## Project at a Glance

| Item | Summary |
| --- | --- |
| Official Documentation | [spatialsnake Documentation](https://spatialsnake-tutorial.readthedocs.io/en/latest/) |
| Tutorial Article | [Core Analysis Tutorial](https://spatialsnake-tutorial.readthedocs.io/en/latest/core_analysis/index.html) |
| Workflow Modes | `single_analysis`, `compare_analysis` |
| Utility Entry Points | `useful_tool`, `produce-file`, `install-packages` |
| Main Analysis Options | `integrate`, `preprocess`, `clustering`, `reclustering`, `annotation_help`, `annotation`, `advance_analysis`, `compare_stage` |
| Supported Input Types | `visium`, `visium_segment`, `visium_HD`, `xenium`, `Merfish`, `stereo_seq` |

## Core Functions

- Standardize raw spatial transcriptomics data into a unified object during `Ingesting`.
- Run `preprocess` for quality control, filtering, normalization, and dimensionality reduction preparation.
- Perform `clustering` and visualization, followed by `annotation_help` and `annotation`.
- Carry out `reclustering` and `reannotation` for clusters of interest.
- Execute `advance_analysis` for downstream analyses and `compare_stage` for cross-sample comparison.
- Use auxiliary utilities for `splitting`, `merge`, and `transform`.

## Available Platforms

### Sequencing-based

- `visium`: 10x Genomics spatial transcriptomics data
- `visium_HD`: high-resolution 10x Genomics spatial transcriptomics data
- `visium_segment`: cell segmentation outputs from 10x Genomics Space Ranger
- `stereo_seq`: BGI Stereo-seq spatial transcriptomics data, including different bin sizes, `cellbin`, and adjusted `cellbin` data types

### Imaging-based

- `xenium`: image-based 10x Genomics Xenium spatial transcriptomics data
- `Merfish`: Vizgen MERFISH spatial transcriptomics data

## Basic Installation

### 1. Create the base conda environment

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n spatialsnake_env python=3.12.11 snakemake-minimal=9.8.1 r-base=4.4.0 -y
conda activate spatialsnake_env
```

### 2. Install the documented core dependencies

```bash
conda install -c conda-forge r-optparse r-tidyverse r-future r-jsonlite r-rcolorbrewer r-patchwork r-cowplot r-pheatmap r-seurat r-remotes r-biocmanager r-presto r-nmf r-circlize
conda install -c bioconda bioconductor-annotationdbi bioconductor-complexheatmap bioconductor-clusterprofiler bioconductor-edger bioconductor-org.hs.eg.db bioconductor-org.mm.eg.db bioconductor-rhdf5 bioconductor-biocneighbors
conda install -c conda-forge bbknn cython
```

### 3. Install `spatialsnake`

#### Option 1. Install from PyPI (recommended)

```bash
pip install spatialsnake
spatialsnake --version
```

#### Option 2. Install from Bioconda

The official Bioconda package provides an alternative Conda-based installation
in a fresh environment:

```bash
conda create -n spatialsnake_env -c conda-forge -c bioconda spatialsnake -y
conda activate spatialsnake_env
spatialsnake --version
spatialsnake install-packages
```

See package details and release files on [Bioconda](https://anaconda.org/bioconda/spatialsnake).

#### Option 3. Install from source code

```bash
git clone https://github.com/zhenghlin/spatialsnake.git
cd spatialsnake
python -m pip install .
python -m pip install ".[extended]"
spatialsnake --version
```

### 4. Optional extended package step

For PyPI or source installs:

```bash
pip install "spatialsnake[extended]"
spatialsnake install-packages
```

For conda installs:

```bash
spatialsnake install-packages --extended
```

With the minimal installation, the documented workflow includes `integrate`, `preprocess`, `clustering`, `reclustering`, `annotation_help`, `annotation`, `reannotation`, and utility operations for merge and split. For conda installs, `spatialsnake install-packages` completes the minimal pip-only core packages; `spatialsnake install-packages --extended` adds downstream Python packages, `pybanksy`, and R/GitHub packages for documented extended components including `compare_stage`, `transform`, `banksy`, and `cellchat`-related workflows. For PyPI installs, keep using `pip install "spatialsnake[extended]"` before `spatialsnake install-packages`.

## Working Directory

Prepare the working directory before running the main workflow:

```text
project_root/
├── data/
├── sample.txt
├── results/
└── <analysis_option>.yaml
```

```bash
mkdir -p project_root/data project_root/results
touch project_root/sample.txt
```

`sample.txt` is the required sample information table for every module in the main workflow. In the working directory, `data/` stores raw input data, `results/` stores analysis outputs generated by the workflow, and `<analysis_option>.yaml` is an optional configuration file.

## Minimal Usage

The command-line interface provides the following documented entry points:

```bash
spatialsnake <command> <INPUT> <TYPE> [--option=<analysis_option>] [options]
spatialsnake useful_tool [--option=<ways>] <INPUT> [options]
spatialsnake produce-file [--option=<analysis_option>]
spatialsnake install-packages [--extended] [--dry-run]
spatialsnake (-h | --help)
spatialsnake --version
```

Main workflow selection:

- `<command>`: choose `single_analysis` or `compare_analysis`
- `<TYPE>`: choose from `visium`, `visium_segment`, `visium_HD`, `xenium`, `Merfish`, and `stereo_seq`
- `--option=<analysis_option>`: choose from `integrate`, `preprocess`, `clustering`, `reclustering`, `annotation_help`, `annotation`, `advance_analysis`, and `compare_stage`

Configuration files can be generated with:

```bash
spatialsnake produce-file --option=<analysis_option>
```

The generated YAML template can then be applied with `--configfile`. Parameters provided directly on the command line take priority over parameters defined in the YAML file.

## Further Reading

- Read the full [official documentation](https://spatialsnake-tutorial.readthedocs.io/en/latest/).
- Start from the example-based [core analysis tutorial](https://spatialsnake-tutorial.readthedocs.io/en/latest/core_analysis/index.html).
- If you encounter problems or would like to suggest extensions, please [open an issue on GitHub](https://github.com/l-zh007/spatialsnake/issues).

## Reference

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., and Nahnsen, S. *Sustainable data analysis with Snakemake*. F1000Research, 10:33, 2021. https://doi.org/10.12688/f1000research.29032.2