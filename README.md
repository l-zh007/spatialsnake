<p align="center">
  <img src="https://raw.githubusercontent.com/l-zh007/spatialsnake/main/spatialsnake-logo.png" alt="Spatialsnake logo" width="700">
</p>

<p align="center"><strong>Reusable spatial transcriptomics workflows for single- and multi-sample analysis</strong></p>

<p align="center">
  <a href="https://spatialsnake-tutorial.readthedocs.io/en/latest/">
    <img src="https://img.shields.io/badge/docs-Read%20the%20Documentation-blue" alt="Documentation">
  </a>
  <a href="https://pypi.org/project/spatialsnake/">
    <img src="https://img.shields.io/pypi/v/spatialsnake" alt="PyPI version">
  </a>
  <a href="https://anaconda.org/bioconda/spatialsnake">
    <img src="https://anaconda.org/bioconda/spatialsnake/badges/version.svg" alt="Bioconda version">
  </a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="MIT license"></a>
</p>

Spatialsnake is a command-line tool for analysing **processed spatial transcriptomics data** from sequencing- and imaging-based platforms. Built on SpatialData and Snakemake, it connects configurable core analysis, targeted reanalysis and downstream methods through reusable SpatialData objects stored in Zarr.

SpatialData provides the Python data model linking expression tables to spatial elements and coordinate systems; Zarr is the on-disk storage format. Spatialsnake coordinates readers, analysis stages and method-specific inputs around these objects. Its ingestion step starts from platform software outputs, such as Space Ranger or SAW results, rather than FASTQ read alignment.

## Key features

- **Core analysis:** ingest supported platform outputs, review QC and dimensionality diagnostics, cluster observations, inspect markers or reference-based estimates, and store reviewed annotations. Compatible samples can share the core workflow while retaining sample and group identities.
- **Focused analysis:** select populations, ROIs or physical-distance neighborhoods; recluster selected observations and return refined labels to the parent object. Prepared subsets can proceed to downstream analyses.
- **Further analyses:** use BANKSY or CellCharter for spatial domains, CellPhoneDB/LIANA/CellChat for communication, and pySCENIC for regulatory programs. Sample-aware comparisons support pseudobulk gene analysis and CellChat comparisons. Supported communication plots can be regenerated from saved inference results.
- **Optional data operations:** split or merge compatible objects, import AnnData, and export AnnData or Seurat objects through the CLI. Saved stage inputs allow selected analyses to be revised without repeating unrelated preparation.

## Supported platforms

The `TYPE` values below are case-sensitive documented spellings. Platform-specific selectors determine the representation loaded for analysis.

| Platform / representation | CLI `TYPE` | Input |
| --- | --- | --- |
| Visium | `visium` | Processed expression matrix and spatial files |
| Visium HD | `visium_HD` | Space Ranger binned outputs; choose `bin_size` in the core sample table |
| Segmented Visium | `visium_segment` | Supported Space Ranger segmentation outputs |
| Stereo-seq | `stereoseq` | Supported SAW outputs; choose bin, cellbin or adjusted cellbin with `input_spec` |
| Xenium | `xenium` | Processed expression/transcript data and available images or cell geometries |
| MERFISH / MERSCOPE | `Merfish` | Supported Vizgen expression, transcript, boundary and image outputs |

Availability of images, transcripts and geometries depends on the platform files and loading settings. Multi-sample analysis is not spatial registration or three-dimensional reconstruction of tissue sections.

## Installation

The [installation guide](docs/installation.md) provides Conda/R dependency commands, followed by the PyPI installation, alternative Bioconda instructions and source installation. **Create and activate that base environment before running:**

```bash
python -m pip install spatialsnake==0.0.6
spatialsnake --version
python -m pip check
```

Core analysis does not require a GPU. Optional modules have additional Python/R dependencies. The [optional installation instructions](docs/installation.md#optional-analyses-pypi-installations) separate CPU and NVIDIA GPU builds and explain the different PyPI and Bioconda installer routes. Installation-time `--torch-backend` and analysis-time `--device` have different purposes.

<details>
<summary>Developers and reproducible source checkouts</summary>

Follow the [source installation instructions](docs/installation.md#3-source-installation-development-or-a-specified-revision). Use the release tag or commit associated with your analysis; the default branch can contain changes absent from a published package. Record local modifications as well as the package version.

</details>

## Multi-sample analysis

`single_analysis` runs the single-sample route; `compare_analysis` provides the multi-sample route and selected group-comparison workflows. For a multi-sample **Visium core analysis**, use one row per sample:

```text
sample_id    input_path            group
Sample_A     data/Sample_A         Control
Sample_B     data/Sample_B         Control
Sample_C     data/Sample_C         Treated
Sample_D     data/Sample_D         Treated
```

Replace the illustrative names and paths with real inputs. Core Visium HD tables additionally require `bin_size`; Stereo-seq tables require `input_spec`, after `group` for the multi-sample route. See the documentation for supported selectors and optional per-sample QC columns. Reference annotation and downstream modules use their own sample-table schemas, so the core table should not be reused without checking the relevant tutorial.

Use `compare_analysis` with the chosen stage, for example `--option=ingestion`. Here, `ingestion` means **data ingestion and preparation of the combined object**, not batch correction. Batch-correction settings are configured separately during preprocessing. Saved combined objects retain sample/group identities; sample-level statistical comparisons require appropriate biological replication.

## Command reference

| Entry point | Purpose |
| --- | --- |
| `single_analysis` / `compare_analysis` | Single-sample / multi-sample workflows |
| `useful_tool` | Splitting, merging and supported data conversion |
| `produce-file` | Generate a YAML template for a selected stage |
| `install-packages` | Complete dependencies for the selected installation route |

```text
spatialsnake <command> <INPUT_FILE> <TYPE> [--option=<analysis_option>] [options]
spatialsnake useful_tool [--option=<operation>] <INPUT>... [options]
spatialsnake produce-file --option=<analysis_option>
spatialsnake install-packages [--extended] [--torch-backend=<auto|cpu|cu128>] [--dry-run]
spatialsnake --help
spatialsnake --version
```

Analysis options include `ingestion`, `preprocess`, `clustering`, `annotation_help`, `annotation`, `reclustering`, `advanced_analysis` and `compare_stage`. Reannotation uses `--option=annotation --anno_algorithm=reannotation`; it is not a separate top-level `--option`. Utility operations are `splitting`, `merge` and `transform`.

The current CLI recommends `ingestion` and `advanced_analysis`; the legacy spellings `integrate` and `advance_analysis` remain supported. For pySCENIC, use `--runpipe=pyscenic`; `pysenic` remains a compatibility alias. Stage names are normalized internally: `produce-file --option=ingestion` creates **`integrate.yaml`**, and `produce-file --option=advanced_analysis` creates **`advance_analysis.yaml`**. Use `spatialsnake --help` to check the options supported by your installed release.

## Documentation and support

- **Core analysis:** [worked single- and multi-sample tutorial](https://spatialsnake-tutorial.readthedocs.io/en/latest/core_analysis/index.html).
- **Focused analysis, utility operations and further analyses:** [full documentation](https://spatialsnake-tutorial.readthedocs.io/en/latest/), with module-specific sample tables and parameters.
- **Problems or feature requests:** [GitHub issues](https://github.com/l-zh007/spatialsnake/issues). Include the software version, installation route, exact command, relevant YAML/sample table and error log. Omit credentials and private data.

## Citation

Report the Spatialsnake version and source revision used in your analysis. Please also cite the underlying SpatialData framework and Snakemake workflow engine:

- Marconato L, Palla G, Yamauchi KA, et al. **SpatialData: an open and universal data framework for spatial omics.** *Nature Methods*. 2025;22:58–62. [doi:10.1038/s41592-024-02212-x](https://doi.org/10.1038/s41592-024-02212-x).
- Mölder F, Jablonski KP, Letcher B, et al. **Sustainable data analysis with Snakemake.** *F1000Research*. 2021;10:33. [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2).

For individual analytical methods, cite the corresponding original publications listed in the module documentation.

## License

Spatialsnake is distributed under the [MIT license](LICENSE).
