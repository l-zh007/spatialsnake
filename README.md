# Snakemake workflow: `<name>`

[![spatialsnake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>)

A Snakemake workflow for `<description>`

- [Snakemake workflow: `<name>`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

Usage:
    spatialsnake single_analysis <INPUT> <TYPE> [--step=STEP] [options]
    spatialsnake compare_analysis <INPUT> <TYPE> [--step=STEP] [options]
    spatialsnake integrate <INPUT> <TYPE> [options]
    spatialsnake segout <INTEGRATED_ZARR> [options]
    spatialsnake --install-packages
    spatialsnake (-h | --help)
    spatialsnake --version

Main Commands:
    single_analysis      Process single spatial transcriptomics dataset (runs all basic steps except advance_analysis by default)
    compare_analysis     Compare multiple spatial transcriptomics datasets
    integrate            Merge multiple ZARR datasets into one integrated dataset
    segout               Split integrated ZARR data back into individual samples

Step Control:
    --step STEP          Run specific analysis step (for single_analysis and compare_analysis)
                         Options: integrate, preprocess, clustering, annotion_help, annotion, advance_analysis
                         Default: run all basic steps (integrate→preprocess→clustering→annotion_help→annotion)

Input Arguments:
    SAMPLE_LIST          Text file containing sample paths (one per line)
    INTEGRATED_ZARR      Integrated ZARR file to split (for segout command)
    annotion_list
    filter_list

Type Arguments:
    visium
    visium_segment
    visium_HD
    xenium
    Merfish
    slide_seq
    Spatial_ATAC

Basic Configuration:
    --configfile FILE    Configuration file in YAML format [default: config.yaml].

Integration Step Options (--step integrate):
    --SAMPLE_LIST FILE        sample list
    --integration-method TEXT   Integration method [default: harmony].
    --cells_boundaries BOOL    xenium key in load in data [default: False].
    --nucleus_boundaries BOOL  xenium key in load in data [default: False].
    --nucleus_labels BOOL      xenium key in load in data [default: False].
    --morphology_mip BOOL      xenium key in load in data [default: False].

Preprocessing Step Options (--step preprocess):
    --annotion_list FILE    for the filter params in different sample
    --min_cells INT         Minimum spots per gene [default: 3].
    --min_genes INT         Minimum genes per spot [default: 200].
    --variable BOOL         Filter the variable spot to analysis [default: False].
    --harmony BOOL          harmony method [default: True].
    --seg_filter BOOL       to seg filter the differnet sample dataset when command compare_anaysis.
    --NEIGHBORS FLOAT       neighbors for pca umap.
Clustering Step Options (--step clustering):
    --resolution FLOAT   Cluster resolution [default: 0.5].
    --cluster_algorithm TEXT Clustering algorithm [default: leiden].
    --tsene BOOL        umap [default:False]
    --MIN_DIST FLOAT    umap_key [default:0.3]
    --SPREAD FLOAT      umap_key [default:1]

Annotation Help Step Options (--step annotion_help):
    --image_slice BOOL        containing marker genes for cell types[default: False].
    --markers_algorithm TEXT       Automatically detect marker genes [default: wilcoxon].
    --shape_type TEXT         Automatically detect marker genes [default: cell_boundaries].
    --image_type TEXT         Automatically detect marker genes [default: hires].
    --spacies TEXT            Automatically detect marker genes [default: human].
    --slice BOOL              params for the image slice to depandent size [default: False].
    --x1 INT
    --x2 INT
    --y1 INT
    --y2 INT
Compare_analyze Step Options (--step compare_analysis)    
    --cell_focus TEXT         celltype you focus to compare in different sample.
    --compare_algorithm TEXT  compare analysys
Annotation Step Options (--step annotion):
    --annotation-file FILE    Annotation file for cell typing (required for annotion step)
    --anno_algorithm TEXT     Annotation method [default: mannul].
    --shape_type TEXT         Automatically detect marker genes [default: cell_boundaries].
    --image_type TEXT         Automatically detect marker genes [default: hires].
    --slice BOOL              params for the image slice to depandent size [default: False].
    --x1 INT
    --x2 INT
    --y1 INT
    --y2 INT
    --max_epochs_reference INT    params for cell2Location model train and test [default: 250].
    --remove_mt BOOL              params for cell2Location model train and test [default: True].
    --N_cells_per_location INT    params for cell2Location model train and test [default: 30].
    --max_epochs_st INT           params for cell2Location model train and test [default: 30000].
    --device TEXT                 cpu or GPU accelerate [default: cuda].
    
Advanced Analysis Step Options (--step advance_analysis):
    --categrory TEXT        Run  which analysis analysis.
    --pyscenic-input FILE   Input file for PySCENIC analysis.
    --pyscenic-db FILE      PySCENIC database directory.
    --pyscenic-feature FILE path for necessary file of pyscenic.
    --pyscenic-tfs FILE     path for necessary file of pyscenic.
    --cell_attr TEXT        cell_id for pyscenic [default: cell_id]
    --workers INT           workers for pyscenic [default: 8].
    --count-data TEXT       gene type for cellPhoneDB [default: hgnc_symbol].
    --threads INT           workers for cellphoneDB [default: 8].
    --output_name           output name for cellPhoneDB [default: Normal].
    
Cell Segmentation Options (applicable to multiple steps):
    --zarr_file FILE        seg with the sample name or region
    

General Options:
    -j INT, --jobs INT   Number of CPU cores [default: 4].
    --output-dir DIR     Output directory [default: results].

Utility Options:
    --install-packages   Install required packages.
    -u, --unlock         Unlock stalled workflow.
    -r, --remove         Remove all output files.
    -d, --dry            Dry run (simulate execution).
    -h, --help           Show this help message.
    --version            Show version.



## Examples:
    Run all basic steps on single sample (default behavior)
    spatialsnake single_analysis samples.txt
    
    Run only preprocessing step
    spatialsnake single_analysis samples.txt --step preprocess
    
    Run annotation with custom annotation file
    spatialsnake single_analysis samples.txt --step annotion --annotation-file celltypes.txt
    
    Run PySCENIC analysis
    spatialsnake single_analysis samples.txt --step advance_analysis --run-pyscenic --pyscenic-db ./databases
    
    Integrate multiple samples
    spatialsnake integrate samples.txt
    
    Split integrated data
    spatialsnake segout integrated_data.zarr


```bash
cd path/to/snakemake-workflow-name
```

Adjust options in the default config file `config/config.yaml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

To run the workflow with **apptainer** / **singularity**, add a link to a container registry in the `Snakefile`, for example `container: "oras://ghcr.io/<user>/<repository>:<version>"` for Github's container registry.
Run the workflow with:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```

## Authors

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO

- Replace `<owner>` and `<repo>` everywhere in the template with the correct user name/organization, and the repository name. The workflow will be automatically added to the [snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/index.html) once it is publicly available on Github.
- Replace `<name>` with the workflow name (can be the same as `<repo>`).
- Replace `<description>` with a description of what the workflow does.
- Update the [deployment](#deployment-options), [authors](#authors) and [references](#references) sections.
- Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment](#deployment-options) options.
- Do not forget to also adjust the configuration-specific `config/README.md` file.
