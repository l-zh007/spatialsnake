#   Spatialsnake  :

[![spatialsnake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>)

A Snakemake workflow for spatial transcriptomics  powered by spatialdata framework

for more detail of the usage of the pipeline 
**please read the [documentation](https://spatialsnake-tutorial.readthedocs.io/en/latest/)**

- [Snakemake workflow: `spatialsnake`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [lzh](#authors)
  - [References](#references)
  - [TODO](#todo)

## ** How to install spatialsnake **


### prepare the environment first.
```bash
## Create conda environment with the environment.yml file in github code page

conda env create -f environment.yml -n spatialsnake_env     ## [or setting your own conda env name]

conda activate spatialsnake_env
```


### install the spatialsnake.

```bash
git clone https://github.com/l-zh007/spatialsnake.git

cd spatialsnake
pip install -e .                       # or [pip install .]  优选-e开发者模式

spatialsnake -h
spatialsnake install-packages          # (Install coordinate R packages)

mkdir project
cd project

start your analysis with file[sample.txt] and spatialdata in [data/*] and output dir [results]

please make sure your spatialdata folder name in [data/] in accordance with sample_name in  [sample.txt]
```


## Examples:
    Run only one step of ["integrate","preprocess","clustering","annotion_help","annotion","compare_analyze","advance_analysis"]
      spatialsnake <sample_channal> sample.txt <data_type> --option=<option_name> [other_params]

    IF you want to run with integrate multiple sample:

    Run compare_analysis
      spatialsnake compare_analysis sample.txt <data_type> --option=<option_name>
      
    Run all basic steps on single sample (default behavior)
      spatialsnake <sample_channal> sample.txt <data_type> --option=all
    
    
    IF you want to run with setting your own params:

      spatialsnake produce-file [--option=<analysis_option>]
      spatialsnake <sample_channal> sample.txt <data_type> --option=<option_name> --config-file <*.yaml>
    
    We produce some useful_tool to help you analysis.
    
    Split integrated data with barcode
      spatialsnake useful_tool integrated_data.zarr --data_barcode B_cell
      
    Split integrated data with image coordinate
      spatialsnake useful_tool integrated_data.zarr --max_x --min_x --max_y --min_y2

    Install coordinate R packages  
      spatialsnake install-packages


## Usage:
    spatialsnake <command> <INPUT> <TYPE> [--option=<analysis_option>] [options]
    spatialsnake useful_tool [--option=<ways>] <INPUT> [options]
    spatialsnake produce-file [--option=<analysis_option>]
    spatialsnake install-packages
    spatialsnake (-h | --help)
    spatialsnake --version

    commands:
      single_analysis      Process single spatial transcriptomics dataset (runs all basic steps except advance_analysis by default)
      compare_analysis     Compare multiple spatial transcriptomics datasets

    analysis option:
        integrate
        preprocess
        clustering
        annotion_help
        annotion
        compare_analyze
        advance_analysis

    Type Arguments:
        visium
        visium_segment
        visium_HD
        xenium
        Merfish
        slide_seq

    INPUT Arguments:
        sample.txt
        annotion.txt
        filter_list
    
    Basic Configuration:
        --configfile <FILE>    Configuration file in YAML format [default: config.yaml].

    Integration Step Options (--option integrate):
        --cells_boundaries <BOOL>    xenium key in load in data [default: False].
        --nucleus_boundaries <BOOL>  xenium key in load in data [default: False].
        --nucleus_labels <BOOL>      xenium key in load in data [default: False].
        --morphology_mip <BOOL>      xenium key in load in data [default: False].

    Preprocessing Step Options (--option preprocess):
        --min_cells <INT>         Minimum spots per gene [default: 3].
        --min_genes <INT>         Minimum genes per spot [default: 200].
        --seg_filter <BOOL>       to seg filter the differnet sample dataset when command compare_anaysis [default: False].
        --filter_list <FILE>      filename of filter [default: False]
        --batch_method <TEXT>     batch method for multiple sample analysis [default: harmony]
        --sketch <BOOL>           whether use sketch method to analysis [default: False]
    
    Clustering Step Options (--option clustering):
        --resolution <FLOAT>        Cluster resolution [default: 0.5].
        --cluster_algorithm <TEXT>  Clustering algorithm [default: leiden].
        --tsene <BOOL>              umap [default:False].
        --n_clusters <INT>          kmeans params of cluster [default: 15].
    
    Annotation Help Step Options (--option annotion_help):
        --markers_algorithm <TEXT>       Automatically detect marker genes [default: wilcoxon].
        --spacies <TEXT>            Automatically detect marker genes [default: human].
    
    Compare_analyze option Options (--option compare_analysis)
        --cell_focus <TEXT>         celltype you focus to compare in different sample[default: None].
        --compare_algorithm <TEXT>  compare analysys [default: DEseq2].
    Annotation option Options (--option annotion):
        --annotation-file <FILE>    Annotation file for cell typing (required for annotion step)
        --anno_algorithm <TEXT>     Annotation method [default: mannul].
        --shape_type <TEXT>         Automatically detect marker genes [default: cell_boundaries].
        --image_type <TEXT>         Automatically detect marker genes [default: hires].
        --device <TEXT>                 cpu or GPU accelerate [default: cuda].

    Advanced Analysis option Options (--option advance_analysis):
        --runpipe <TEXT>        Run  which analysis analysis.[default: advance_analysis]
        --senic_input <DIR>   Input file for PySCENIC analysis.[default: sample.zarr]
        --motifs_input <FILE>      PySCENIC database directory.[default: motifs-v9-nr.hgnc-m0.001-o0.0.tbl]
        --feather_input <FILE> path for necessary file of pyscenic.[default: hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather]
        --tfs_input <FILE>     path for necessary file of pyscenic.[default: hs_hgnc_tfs.txt]
        --count-data <TEXT>       gene type for cellPhoneDB [default: hgnc_symbol].
        --threads <INT>           workers for cellphoneDB [default: 8].
        --output_name <TEXT>      output name for cellPhoneDB [default: Normal].
    useful_tool params Option:
        --output_zarr_path <FILE> output dir for splitted file [default: results]
        --split_by <TEXT>          slice out with the barcode in table[anndata] .obs [default: clusters]
        --max_x   <FLOAT>         coordinate of image boundaries [default: 0]
        --min_x   <FLOAT>         coordinate of image boundaries [default: 2000]
        --max_y   <FLOAT>         coordinate of image boundaries [default: 2000]
        --min_y   <FLOAT>         coordinate of image boundaries [default: 0]
    
    General Options:
        -j <INT>, --jobs <INT>   Number of CPU cores [default: 32].
        --results_folder <DIR>     Output directory [default: results].

    Utility Options:
        --install-packages   Install required packages.
        -u, --unlock         Unlock stalled workflow.
        -r, --remove         Remove all output files.
        -d, --dry            Dry run (simulate execution).
        -h, --help           Show this help message.
        --version            Show version.


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
