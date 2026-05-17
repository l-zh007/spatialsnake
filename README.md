<table border="0" cellspacing="0" cellpadding="0">
  <tbody>
    <tr>
      <td rowspan="3" valign="top" style="padding-right: 20px;">
        <h1>spatialsnake</h1>
        <p>A Snakemake workflow for spatial transcriptomics powered by the <code>spatialdata</code> framework.</p>
      </td>
      <td>
        <a href="https://spatialsnake-tutorial.readthedocs.io/en/latest/">
          <img src="https://img.shields.io/badge/docs-Read%20the%20Documentation-blue" alt="Documentation">
        </a>
      </td>
    </tr>
    <tr>
      <td>
        <a href="https://pypi.org/project/spatialsnake/">
          <img src="https://img.shields.io/pypi/v/spatialsnake" alt="PyPI">
        </a>
      </td>
    </tr>
    <tr>
      <td>
        <a href="https://snakemake.github.io">
          <img src="https://img.shields.io/badge/snakemake-%E2%89%A58.0.0-brightgreen.svg" alt="Snakemake">
        </a>
      </td>
    </tr>
  </tbody>
</table>

`spatialsnake` provides a command-line workflow for single-dataset analysis, multiple-dataset comparison, stepwise execution of core analysis stages, and auxiliary utilities for data splitting, merging, and format transformation.

## Project at a Glance

| Item | Link / Summary |
| --- | --- |
| Documentation | [spatialsnake Documentation](https://spatialsnake-tutorial.readthedocs.io/en/latest/) |
| Tutorial Article | [spatialsnake Tutorial](https://spatialsnake-tutorial.readthedocs.io/en/latest/) |
| Repository | [GitHub Project Page](https://github.com/l-zh007/spatialsnake) |
| Core Commands | `single_analysis`, `compare_analysis`, `useful_tool`, `produce-file`, `install-packages` |
| Core Analysis Options | `integrate`, `preprocess`, `clustering`, `reclustering`, `annotation_help`, `annotation`, `compare_stage`, `advance_analysis` |

## Quick Installation

```bash
conda env create -f environment.yml -n spatialsnake_env
conda activate spatialsnake_env
git clone https://github.com/l-zh007/spatialsnake.git
cd spatialsnake
pip install -e .
spatialsnake -h
spatialsnake install-packages
```

## Core Functions

- Process a single spatial transcriptomics dataset with `single_analysis`.
- Compare multiple spatial transcriptomics datasets with `compare_analysis`.
- Run analysis stages independently, including integration, preprocessing, clustering, annotation support, annotation, comparative analysis, reclustering, and advanced analysis.
- Use auxiliary tools for `splitting`, `merge`, and `transform`.
- Work with the following input types: `visium`, `visium_segment`, `visium_HD`, `xenium`, `Merfish`, `slide_seq`, and `stereoseq`.

## Working Directory

### Prepare a working directory

```bash
mkdir project
cd project
```

Start the analysis with `sample.txt`, spatialdata stored in `data/*`, and the output directory `results`.

Please ensure that each spatialdata folder name under `data/` is consistent with the corresponding `sample_name` recorded in `sample.txt`.

## Minimal Usage

### Run one analysis step

```bash
spatialsnake <command> sample.txt <TYPE> --option=<analysis_option> [options]
```

### Run multiple-sample comparison

```bash
spatialsnake compare_analysis sample.txt <TYPE> --option=<analysis_option>
```

### Run all basic steps for a single dataset

```bash
spatialsnake single_analysis sample.txt <TYPE> --option=all
```

### Generate a configuration file

```bash
spatialsnake produce-file [--option=<analysis_option>]
spatialsnake <command> sample.txt <TYPE> --option=<analysis_option> --configfile <FILE>
```

### Run utility tools

```bash
spatialsnake useful_tool [--option=<ways>] <INPUT>... [options]
```

## Command Summary

### Workflow commands

- `single_analysis`: process a single spatial transcriptomics dataset; by default, all basic steps except `advance_analysis` are executed.
- `compare_analysis`: compare multiple spatial transcriptomics datasets.

### Analysis options

- `integrate`
- `preprocess`
- `clustering`
- `reclustering`
- `annotation_help`
- `annotation`
- `compare_stage`
- `advance_analysis`

### Utility commands

- `produce-file`: generate a configuration file for a selected analysis option.
- `install-packages`: install required packages.
- `useful_tool --option=splitting`: split integrated data.
- `useful_tool --option=merge`: merge data.
- `useful_tool --option=transform`: transform data formats.

## Selected Parameters

### Basic configuration

- `--configfile <FILE>`: configuration file in YAML format. Default: `config.yaml`.
- `-j, --jobs <INT>`: number of CPU cores. Default: `16`.
- `--results_folder <DIR>`: output directory. Default: `results`.

### Integration

- `--cells_boundaries <BOOL>`
- `--nucleus_boundaries <BOOL>`
- `--nucleus_labels <BOOL>`
- `--morphology_mip <BOOL>`
- `--bin_size <INT>`
- `--merscope_z_layers <TEXT>`
- `--merscope_region_name <TEXT>`

### Preprocessing

- `--min_cells <INT>`
- `--min_genes <INT>`
- `--seg_filter <BOOL>`
- `--filter_list <FILE>`
- `--batch_method <TEXT>`
- `--sketch <BOOL>`
- `--mt_threshold <FLOAT>`

### Clustering and annotation

- `--resolution <FLOAT>`
- `--cluster_algorithm <TEXT>`
- `--n_clusters <INT>`
- `--pcs <INT>`
- `--markers_algorithm <TEXT>`
- `--spacies <TEXT>`
- `--annotation-file <FILE>`
- `--anno_algorithm <TEXT>`

### Advanced analysis and utilities

- `--runpipe <TEXT>`
- `--senic_input <DIR>`
- `--motifs_input <FILE>`
- `--feather_input <FILE>`
- `--tfs_input <FILE>`
- `--count-data <TEXT>`
- `--threads <INT>`
- `--output_name <TEXT>`
- `--split_by=<TEXT>`
- `--output_dir=<TEXT>`

For detailed parameter descriptions, please consult the [official documentation](https://spatialsnake-tutorial.readthedocs.io/en/latest/).

## Reference

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., and Nahnsen, S. *Sustainable data analysis with Snakemake*. F1000Research, 10:33, 2021. https://doi.org/10.12688/f1000research.29032.2
