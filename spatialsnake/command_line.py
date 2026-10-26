#!/usr/bin/env python
'''
spatialsnake
@author: Zhenghao Lin
@email: l-zh007@users.noreply.github.com
'''
import datetime
import importlib.metadata
import importlib.util
import json
import logging
import os
import random
import shlex
import shutil
import subprocess
import sys
import tempfile
import timeit
import warnings
from typing import Any, Dict, List, Optional, Tuple

try:
    import resource
except ImportError:
    resource = None

import spatialsnake
import yaml
from docopt import docopt
from spatialsnake.workflow.function.config_validator import (
    ConfigValidationError,
    OPTION_ALIASES,
    RUNPIPE_ALIASES,
    validate_merge_config,
    validate_splitting_config,
    validate_transform_config,
    validate_workflow_config,
)
from yaml.loader import SafeLoader

# Filter warnings
warnings.filterwarnings("ignore")

# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s - %(message)s',
    stream=sys.stdout
)
logger = logging.getLogger(__name__)
LOG_WIDTH = 80


def cli_option_provided(option_name: str) -> bool:
    """Return True only when the user explicitly typed an option."""
    return any(arg == option_name or arg.startswith(f"{option_name}=") for arg in sys.argv)


def normalize_cli_equals_tokens(argv: List[str]) -> List[str]:
    """Accept ``--option = value`` in addition to standard docopt syntax."""
    normalized = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if (
            token.startswith("--")
            and index + 2 < len(argv)
            and argv[index + 1] == "="
        ):
            normalized.append(f"{token}={argv[index + 2]}")
            index += 3
            continue
        normalized.append(token)
        index += 1
    return normalized


CLI_KEY_ALIASES = {
    "annotation-file": "annotation_list",
    "count-data": "counts_data",
    "seg_filter": "filter_list",
    # ``threads`` is the canonical workflow-wide concurrency setting.  Keep
    # the historical module-specific flags as CLI aliases so existing command
    # lines continue to work without allowing individual scripts to silently
    # oversubscribe the Snakemake allocation.
    "workers": "threads",
    "senic_workers": "threads",
    "cellchat_workers": "threads",
    "max_cores": "threads",
}


def config_key_from_cli(option_name: str) -> str:
    """Map a docopt option name to the config key used by yaml/Snakemake."""
    clean_key = option_name.lstrip("--")
    return CLI_KEY_ALIASES.get(clean_key, clean_key)

# Constants
SPATIALSNAKE_PATH = os.path.dirname(spatialsnake.__file__)
VALID_OPTIONS = [
    "integrate", "preprocess", "clustering", "annotation_help", "annotation",
    "compare_stage", "advance_analysis", "reclustering", "splitting", "merge", "transform"
]

__author__ = 'lzh'
__version__ = '0.0.7'
__logo__ = """

  ╭─── SpatialSnake · v0.0.7 ───╮
  │                             │
  │    ●───●───●───●───●───●    │
  │    │ ╲ │ ╱ │ ╲ │ ╱ │ ╲ │    │
  │    ●───●───●───●───●───●    │
  │                             │
  │   Spatial Transcriptomics   │
  │      Analysis Pipeline      │
  ╰─────────────────────────────╯
"""

__licence__ = """
MIT License
Copyright (c) 2025
...
"""

__doc__ = f"""Main spatialsnake executable, version: {__version__}
{__logo__} 

Usage:
    spatialsnake useful_tool [--option=<ways>] <INPUT>... [options]
    spatialsnake <command> <INPUT_FILE> <TYPE> [--option=<analysis_option>] [options]
    spatialsnake produce-file [--option=<analysis_option>]
    spatialsnake install-packages [--extended] [--torch-backend=<backend>] [--dry-run]
    spatialsnake (-h | --help)
    spatialsnake --version

commands:
    single_analysis      Process single spatial transcriptomics dataset (runs all basic steps except advanced_analysis by default)
    compare_analysis     Run multi-sample spatial transcriptomics analysis

analysis option:
    ingestion           Data ingestion and QC (legacy alias: integrate)
    preprocess
    clustering
    reclustering
    annotation_help
    annotation
    compare_stage
    advanced_analysis   Advanced analysis (legacy alias: advance_analysis)

Type Arguments:
    visium
    visium_segment
    visium_HD
    xenium
    Merfish
    stereoseq

INPUT Arguments:
    sample.txt
    annotation.txt
    filter_list
    
Basic Configuration:
    --configfile <FILE>    User YAML configuration file. If omitted, Spatialsnake uses the matching workflow/envs template.
                           YAML values are authoritative defaults; command-line values override only when explicitly provided.

Data Ingestion and QC Options (--option ingestion; integrate remains supported):
    --cells_boundaries <BOOL>    xenium key in load in data (yaml default: integrate.yaml).
    --nucleus_boundaries <BOOL>  xenium key in load in data (yaml default: integrate.yaml).
    --nucleus_labels <BOOL>      xenium key in load in data (yaml default: integrate.yaml).
    --morphology_mip <BOOL>      xenium key in load in data (yaml default: integrate.yaml).
    --merscope_z_layers <TEXT>   optional z layers for `spatialdata_io.merscope`, e.g. `0` or `0,1,2`.
    --merscope_region_name <TEXT> optional region name for `spatialdata_io.merscope`.
    --merscope_transcripts <BOOL> load transcripts in `spatialdata_io.merscope` (yaml default: integrate.yaml).
    --merscope_cells_boundaries <BOOL> load cell boundaries in `spatialdata_io.merscope` (yaml default: integrate.yaml).
    --merscope_cells_table <BOOL> load cells table in `spatialdata_io.merscope` (yaml default: integrate.yaml).
    --merscope_mosaic_images <BOOL> load mosaic images in `spatialdata_io.merscope` (yaml default: integrate.yaml).

Preprocessing Step Options (--option preprocess):
    --min_cells <INT>         Minimum spots/cells per gene (yaml default: preprocess.yaml).
    --min_counts <INT>        Minimum total counts per spot/cell (yaml default: preprocess.yaml).
    --seg_filter <BOOL>       Alias of --filter_list for YAML per-sample QC overrides (yaml default: preprocess.yaml).
    --filter_list <BOOL>      Enable YAML per-sample QC overrides; sample.txt filter columns are detected automatically.
    --batch_method <TEXT>     Batch correction method for multiple sample analysis (yaml default: preprocess.yaml).
    --sketch <BOOL>           Whether use sketch method to analysis (yaml default: preprocess.yaml).
    --variable <BOOL>         Use highly variable genes for PCA (yaml default: preprocess.yaml).
    --mt_threshold <FLOAT>    Mitochondrial percent threshold for filtering (yaml default: preprocess.yaml).
    --n_top_genes <INT>       Number of highly variable genes (yaml default: preprocess.yaml).
    --n_comps <INT>           Number of PCA components (yaml default: preprocess.yaml or clustering.yaml).
    --NEIGHBORS <INT>         Number of neighbors for graph construction (yaml default: preprocess.yaml or clustering.yaml).
    --sample_rate <FLOAT>     Sketch/downsample rate in (0, 1] (yaml default: preprocess.yaml).
    
Clustering Step Options (--option clustering):
    --resolution <FLOAT>        Cluster resolution (yaml default: clustering.yaml).
    --cluster_algorithm <TEXT>  Clustering algorithm (yaml default: clustering.yaml).
    --tsne <BOOL>               Whether to compute t-SNE in addition to UMAP (yaml default: clustering.yaml).
    --n_clusters <INT>          K-means cluster number (yaml default: clustering.yaml).
    --pcs <INT>                 PCA dimensions used for neighbors (yaml default: clustering.yaml).
    --MIN_DIST <FLOAT>          UMAP min_dist parameter (yaml default: clustering.yaml).
    --SPREAD <FLOAT>            UMAP spread parameter (yaml default: clustering.yaml).

Reclustering Step Options (--option reclustering):
    --recluster_resolution <FLOAT>      Leiden resolution (yaml default: reclustering.yaml).
    --recluster_n_top_genes <INT>       0 uses all genes without HVG selection (default); >=2 selects subset HVGs.
    --recluster_neighbors <INT>         Neighbors for graph construction (yaml default: reclustering.yaml).
    --recluster_n_pcs <INT>             Number of PCs for neighbors (yaml default: reclustering.yaml).
    --recluster_marker_method <TEXT>    Marker test method (yaml default: reclustering.yaml).
    --recluster_min_pct <FLOAT>         Min fraction for marker filtering (yaml default: reclustering.yaml).
    --recluster_logfc_threshold <FLOAT> Min log2FC for marker filtering (yaml default: reclustering.yaml).
    
Annotation Support Options (--option annotation_help):
    --markers_algorithm <TEXT>       Marker gene test method (yaml default: annotation_help.yaml).
    --species <TEXT>                 Species for enrichment or CellPhoneDB: human or mouse (yaml default: module config).

Compare_stage option Options (--option compare_stage)
    sample.txt                 Three whitespace-separated columns for compare_gene: sample_id, input_path, and group.
    --cell_focus <TEXT>         Cell type to compare; accepts one value, comma-separated values, or all (yaml default: compare_stage.yaml).
    --compare_algorithm <TEXT>  Formal pseudobulk DEG algorithm: DESeq2 or edgeR (yaml default: compare_stage.yaml).
    --compare_input_zarr <DIR>  Optional annotated merged zarr override (default: results/merge_data/annotation/concatenated_sdata.zarr).
    --compare_contrasts <TEXT>  Optional comparison:reference pairs; omitted for two groups, required for three or more groups.
    --compare_celltype_col <TEXT> Cell type column in zarr obs (yaml default: compare_stage.yaml).
    --compare_sample_col <TEXT> Sample/biological replicate column in zarr obs (yaml default: compare_stage.yaml).
    --compare_condition_col <TEXT> Optional obs group column checked for stale labels; sample.txt remains authoritative.
    --count_layer <TEXT>        Preferred raw-count layer; standard raw-count sources are validated as fallback.
    --min_replicates <INT>      Minimum biological replicates per condition (yaml default: compare_stage.yaml).
    --min_cells_per_sample <INT> Minimum cells/spots per sample in each cell type (yaml default: compare_stage.yaml).
    --min_total_counts_per_gene <INT> Minimum total pseudobulk counts per gene (yaml default: compare_stage.yaml).
    --de_top_n <INT>            Number of top genes highlighted in DEG figures (yaml default: compare_stage.yaml).
    --cut_off_pvalue <FLOAT>    Adjusted p-value cutoff for DEG figures and enrichment (yaml default: compare_stage.yaml).
    --cut_off_logFC <FLOAT>     Absolute log2FC cutoff for DEG figures and split tables (yaml default: compare_stage.yaml).
    --compare_gene_id_type <TEXT> Gene IDs used for enrichment: auto, SYMBOL, ENSEMBL, or ENTREZID (yaml default: compare_stage.yaml).
    --compare_gene_symbol_col <TEXT> Optional zarr var column containing gene symbols (yaml default: compare_stage.yaml).
    --compare_go_ontology <TEXT> GO ontology: BP, CC, MF, or ALL (yaml default: compare_stage.yaml).
    --compare_enrichment_top_n <INT> Maximum enriched terms displayed per panel (yaml default: compare_stage.yaml).
    --cellchat_compare_output_dir <DIR> Output directory for two-condition CellChat comparison.
    --cellchat_compare_focus_cells <TEXT> Compare all directions among comma-separated cell groups.
    --cellchat_compare_cell_pairs <TEXT> Exact CellChat directions, e.g. A|B,A<->C; overrides other cell scopes.
    --cellchat_compare_source_cells <TEXT> Optional comma-separated sender groups.
    --cellchat_compare_target_cells <TEXT> Optional comma-separated receiver groups.
    --cellchat_compare_pathways <TEXT> Optional pathways; empty selects differential pathways automatically.
    --cellchat_compare_lr_pairs <TEXT> Optional ligand|receptor pairs; empty selects differential LR automatically.
    --cellchat_compare_top_cell_pairs <INT> Number of directed cell pairs selected automatically.
    --cellchat_compare_top_pathways <INT> Number of pathways selected automatically within the cell scope.
    --cellchat_compare_top_lr <INT> Maximum unique LR interactions shown in focused comparison plots.
    --cellchat_compare_plot_advanced <BOOL> Generate selected-pathway official CellChat plots.
Annotation option Options (--option annotation):
    --annotation-file <FILE>    Annotation file for cell typing (required for annotation step)
    --anno_algorithm=<METHOD>   Annotation method: manual, reannotation, cell2Location, or RCTD. Use manual for manual annotation.
    --shape_type <TEXT>         Shape element used for spatial visualization (yaml default: annotation.yaml or annotation_help.yaml).
    --image_type <TEXT>         Image element used for spatial visualization (yaml default: annotation.yaml or annotation_help.yaml).
    --vis_mode <TEXT>           Spatial visualization mode: auto/point/shape (yaml default: annotation.yaml or annotation_help.yaml).
    --point_size <FLOAT>        Point render size for point-based visualization (yaml default: annotation.yaml or annotation_help.yaml).
    --image_slice <BOOL>        Crop spatial visualization to x/y bounds (yaml default: annotation.yaml or annotation_help.yaml).
    --x1 <FLOAT>                Crop minimum x coordinate (yaml default: annotation.yaml or annotation_help.yaml).
    --x2 <FLOAT>                Crop maximum x coordinate (yaml default: annotation.yaml or annotation_help.yaml).
    --y1 <FLOAT>                Crop minimum y coordinate (yaml default: annotation.yaml or annotation_help.yaml).
    --y2 <FLOAT>                Crop maximum y coordinate (yaml default: annotation.yaml or annotation_help.yaml).
    --device <TEXT>             Compute device: auto, cpu, or cuda (yaml default: selected module).
    --max_cores <INT>           Legacy alias of --threads for annotation tools.
    --RCTD_mode <TEXT>             RCTD mode: full or doublet (yaml default: annotation.yaml).
Advanced Analysis Options (--option advanced_analysis; advance_analysis remains supported):
    --runpipe <TEXT>          Select module: cellPhoneDB, pyscenic (legacy: pysenic), liana, cellcharter, banksy, or cellchat.
    --senic_input <DIR>       Input file for PySCENIC analysis.
    --motifs_input <FILE>     PySCENIC database directory.
    --feather_input <FILE>    Comma-separated cisTarget ranking feather database path(s) for pySCENIC.
    --tfs_input <FILE>        path for necessary file of pyscenic.
    --pyscenic_top_regulons <INT> Number of pySCENIC regulons shown in dotplot/violin plots.
    --pyscenic_min_regulon_genes <INT> Minimum target genes retained per pySCENIC regulon.
    --count-data <TEXT>       gene type for CellPhoneDB; alias of --counts_data (yaml default: advance_analysis.yaml).
    --counts_data <TEXT>      gene type for CellPhoneDB (yaml default: advance_analysis.yaml).
    --threads <INT>           Threads allocated to each workflow rule. Explicit CLI values override YAML defaults.
    --output_name <TEXT>      output name for modules that support a named result (yaml default: selected module).
    --workers <INT>           Legacy alias of --threads.
    --senic_workers <INT>     Legacy alias of --threads for pySCENIC.
    --threshold <FLOAT>       CellPhoneDB expression proportion threshold (yaml default: advance_analysis.yaml).
    --pvalue <FLOAT>          CellPhoneDB p-value threshold (yaml default: advance_analysis.yaml).
    --iterations <INT>        CellPhoneDB permutation iterations (yaml default: advance_analysis.yaml).
    --celltype_col <TEXT>     Cell type column in zarr/table/obs (yaml default: advance_analysis.yaml or annotation.yaml).
    --niche_col <TEXT>        niche column in the zarr/table/obs (yaml default: advance_analysis.yaml).
    --cell_pairs <TEXT>       Optional CellPhoneDB cell pairs for plotting, e.g. A|B,A<->C (yaml default: advance_analysis.yaml).
    --cpdb_pathway <TEXT>     Optional CellPhoneDB pathway/classification used for focused plotting.
    --interaction_pairs <TEXT> Comma-separated CellPhoneDB ligand-receptor pairs for focused plotting.
    --cpdb_genes <TEXT>       Comma-separated genes for the CellPhoneDB expression summary.
    --liana_method <TEXT>     LIANA method (yaml default: advance_analysis.yaml).
    --liana_resource_name <TEXT> LIANA ligand-receptor resource (yaml default: advance_analysis.yaml).
    --liana_expr_prop <FLOAT> LIANA expression proportion filter (yaml default: advance_analysis.yaml).
    --liana_min_cells <INT>   LIANA minimum cells per cell type (yaml default: advance_analysis.yaml).
    --liana_top_n <INT>       Number of LIANA interactions highlighted in plots (yaml default: advance_analysis.yaml).
    --liana_pvalue <FLOAT>    LIANA p-value cutoff for p-value-based plotting (yaml default: advance_analysis.yaml).
    --liana_source_celltypes <TEXT> Optional comma-separated LIANA sender cell types for focused plots.
    --liana_target_celltypes <TEXT> Optional comma-separated LIANA receiver cell types for focused plots.
    --liana_cell_pairs <TEXT> Optional LIANA directed or bidirectional cell pairs, e.g. A|B,A<->C.
    --liana_pairs <TEXT>      Optional comma-separated LIANA ligand|receptor pairs for focused plots.
    --cellchat_species <TEXT> CellChat species (yaml default: advance_analysis.yaml).
    --cellchat_min_cells <INT> CellChat minimum cells per group (yaml default: advance_analysis.yaml).
    --cellchat_future_max_size_gb <FLOAT> Maximum future global size in GB for CellChat; 0 means unlimited.
    --cellchat_workers <INT>  Legacy alias of --threads for CellChat.
    --cellchat_is_single_cell <BOOL> Run CellChat in single-cell mode without spatial distance constraints.
    --cellchat_trim <FLOAT>   CellChat truncated mean trim value (yaml default: advance_analysis.yaml).
    --cellchat_interaction_length <INT> CellChat spatial interaction length (yaml default: advance_analysis.yaml).
    --cellchat_spot_size <FLOAT> CellChat spatial spot/cell-size proxy (yaml default: advance_analysis.yaml).
    --cellchat_db_subset <TEXT> CellChat database subset: all_interactions, secreted_only, or secreted_ecm.
    --cellchat_focus_cells <TEXT> Recommended comma-separated cell groups for focused all-direction CellChat plots.
    --cellchat_cell_pairs <TEXT> Optional exact CellChat pairs, e.g. A|B,A<->C; overrides other cell-scope parameters.
    --cellchat_source_cells <TEXT> Optional comma-separated CellChat sender groups for directional focused plots.
    --cellchat_target_cells <TEXT> Optional comma-separated CellChat receiver groups for directional focused plots.
    --cellchat_pathways <TEXT> Optional comma-separated pathways; empty selects top pathways within the cell scope.
    --cellchat_lr_pairs <TEXT> Optional comma-separated ligand|receptor pairs; empty selects top LR within the cell scope.
    --cellchat_top_cell_pairs <INT> Number of strongest directed cell pairs selected automatically.
    --cellchat_top_pathways <INT> Number of top CellChat pathways selected automatically within the cell scope.
    --cellchat_bubble_top_lr <INT> Maximum unique LR interactions shown in focused bubble plots; 0 shows all.
    --cellchat_plot_advanced <BOOL> Generate selected-pathway aggregate, contribution, expression, and signaling-role plots.
    --cellchat_pair_lr_use <TEXT> Legacy optional ligand-receptor pair for the CellChat spatial LR plot.
    --condition_col <TEXT>    Condition column for CellCharter (yaml default: advance_analysis.yaml).
    --sample_col <TEXT>       Observation column containing sample identities (yaml default: selected module).
    --cellcharter_col <TEXT>  Output column for CellCharter domains (yaml default: advance_analysis.yaml).
    --k_geom <INT>            BANKSY spatial neighbor parameter (yaml default: advance_analysis.yaml or clustering.yaml).
    --max_m <FLOAT>           BANKSY harmonic order parameter (yaml default: advance_analysis.yaml or clustering.yaml).
    --banksy_n_comps <TEXT>   BANKSY PCA dimensions, e.g. 20 or 20,50 (yaml default: advance_analysis.yaml).
    --banksy_resolution <TEXT> BANKSY Leiden resolutions, e.g. 0.5 or 0.4,0.8 (yaml default: advance_analysis.yaml).
    --banksy_num_nn <INT>     BANKSY expression-neighbor count for Leiden graph (yaml default: advance_analysis.yaml).
    --banksy_max_features <INT> Maximum genes used to build dense BANKSY matrices; 0 uses all genes.
    --banksy_feature_col <TEXT> Optional var column used before variance-based feature selection.
    --banksy_add_umap <BOOL>  Compute UMAP for BANKSY matrices (yaml default: advance_analysis.yaml).
    --banksy_plot_full <BOOL> Generate pyBANKSY full result figures (yaml default: advance_analysis.yaml).
    --banksy_run_nonspatial <BOOL> Run non-spatial baseline clustering (yaml default: advance_analysis.yaml).
    --banksy_plot_celltype_enrichment <BOOL> Plot celltype enrichment when obs['celltype'] exists.
    --banksy_plot_max_points <INT> Maximum points shown in lightweight BANKSY spatial plot.
    --banksy_sample_col <TEXT> Sample/slice column used to avoid cross-sample BANKSY neighbors.
    --banksy_selected_lambda <TEXT> Optional final BANKSY lambda written to spatial_cluster.
    --banksy_selected_resolution <TEXT> Optional final BANKSY resolution written to spatial_cluster.
    --banksy_seed <INT>       Random seed for BANKSY Leiden clustering and plotting downsampling.

useful_tool splitting Option:
    --output_dir=<TEXT>       final useful-tool output directory; defaults use a tool-specific folder under results/useful_results.
    --split_by=<TEXT>         slice out with the barcode in table[anndata] .obs (yaml default: splitting.yaml).
    --barcodes=<TEXT>         selected labels for split_by filter: comma combines labels, pipe creates separate outputs (yaml default: splitting.yaml).
    --roi_csv=<TEXT>          csv file or directory for ROI splitting (yaml default: splitting.yaml).
    --shape_elements=<TEXT>   slice out with the shape (yaml default: splitting.yaml).
    --table_key=<TEXT>        SpatialData table used by splitting or merge; empty requires an unambiguous table.
    --subset_mode=<TEXT>      Subset only the table or also associated spatial elements: table or spatial.
    --coordinate_system=<TEXT> Coordinate system for ROI polygons or spatial cropping; empty enables unambiguous inference.
    --roi_label_col=<TEXT>    Optional ROI/category column in Loupe or Xenium CSV files.
    --roi_region=<TEXT>       Optional region used to disambiguate repeated instance IDs in integrated objects.
    --annotation_format=<TEXT> Ordinary splitting annotation CSV: auto, loupe, xenium, both, or none.
    --annotation_cols=<TEXT>  Ordinary splitting obs columns exported as software annotations.
    --distance=<FLOAT>        Outer distance from ROI/target cells in micrometres; enables distance splitting.
    --distance_min=<FLOAT>    Inner distance in micrometres; selected neighbors satisfy min <= d < max.
    --neighbor_types=<TEXT>   Optional comma-separated labels allowed among neighboring observations.
    --neighbor_col=<TEXT>     obs column containing neighbor_types labels.
    --max_x=<FLOAT>           coordinate of image boundaries (yaml default: splitting.yaml).
    --min_x=<FLOAT>           coordinate of image boundaries (yaml default: splitting.yaml).
    --max_y=<FLOAT>           coordinate of image boundaries (yaml default: splitting.yaml).
    --min_y=<FLOAT>           coordinate of image boundaries (yaml default: splitting.yaml).

useful_tool merge Option:
    --merge_by=<TEXT>          merge independent samples or overlay annotations: sample or annotation.
    --sample_ids=<TEXT>        optional comma-separated source IDs in input order.
    --feature_join=<TEXT>      gene join for sample mode: inner or outer.
    --label_cols=<TEXT>        optional comma-separated obs label columns to transform.
    --label_policy=<TEXT>      preserve, prefix, or offset labels during sample merge.
    --annotation_csv=<TEXT>    csv path, directory, or comma-separated csv paths for annotation mode.
    --annotation_col=<TEXT>    annotation column in CSV or subset Zarr; auto enables unambiguous inference.
    --csv_cell_col=<TEXT>      cell ID column in annotation CSV; auto enables inference.
    --csv_region_col=<TEXT>    region/sample column in annotation CSV; required for repeated IDs.
    --input_cell_col=<TEXT>    cell ID column in parent table; auto uses instance_key.
    --input_region_col=<TEXT>  region column in parent table; auto uses region_key.
    --target_col=<TEXT>        parent obs column updated by annotation mode.
    --fallback_col=<TEXT>      parent obs column used to initialize a new target column.
    --conflict_policy=<TEXT>   annotation conflict policy: error, first, or last.
    --existing_policy=<TEXT>   overwrite_matched, fill_missing, or error for an existing target.
    --min_match_rate=<FLOAT>   minimum source-to-parent annotation ID match rate.

useful_tool transform Option:
    --save_image=<BOOL>        preserve downsampled legacy images where supported.
    --transform_from=<TEXT>    source format: zarr or h5ad.
    --transform_to=<TEXT>      target format: zarr, h5ad, or seurat.
    --type=<TEXT>              data type for Seurat export: auto, st, or sc.
    --seurat_matrix=<TEXT>     Seurat matrix source: auto, raw, or X.
    --memory_limit_gb=<FLOAT>  memory safety limit; 0 uses 80 percent of available memory.
    --keep_intermediate=<BOOL> retain the H5AD intermediate for zarr to Seurat.

General Options:
    -j <INT>, --jobs <INT>   Number of CPU cores [default: 16].
    --results_folder <DIR>     Output directory (yaml default: selected workflow template).

Utility Options:
    --install-packages   Install required packages.
    --extended           Install full downstream packages for conda/Bioconda environments.
    --torch-backend <TEXT> Bioconda PyTorch build: auto reuses preinstalled Torch; cpu/cu128 explicitly install that build [default: auto].
    --dry-run            Preview install-packages actions without installing.
    -u, --unlock         Unlock stalled workflow.
    -r, --remove         Remove all output files.
    -d, --dry            Dry run (simulate execution).
    -h, --help           Show this help message.
    --version            Show version.

"""


class BaseCommandLine:
    """Base class for handling command line execution and logging."""

    def __init__(self, arguments: Dict[str, Any]):
        self.arguments = arguments
        self.runid = "".join(random.choices("abcdefghisz", k=3) + random.choices("123456789", k=5))
        self.config: List[str] = []
        self.parameters: Dict[str, Any] = {}
        self.log_enabled = True
        self.cmd_str = ""
        self.configfile_loaded = False
        self.configfile_path = ""
        self.logname = ""
        self.log_header_written = False
        self.execution_status = "NOT_STARTED"
        self.return_code = None
        self.error_message = ""
        self.started_at = None
        self.finished_at = None
        self.resource_before = None
        self.resource_after = None

    def __str__(self):
        return self.cmd_str

    def __repr__(self):
        return self.cmd_str

    def add_config_argument(self):
        """Append config arguments to command string."""
        raise NotImplementedError

    def prepare_arguments(self):
        """Parse arguments and prepare command string."""
        raise NotImplementedError

    def _log_label(self):
        label = self.arguments.get("--option")
        if not label and self.arguments.get("useful_tool"):
            label = f"useful_tool_{self.arguments.get('--option', 'run')}"
        if not label:
            label = self.arguments.get("<command>", "run")
        return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(label or "run"))

    def _ensure_logname(self):
        if not self.logname:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            log_dir = "log"
            self.logname = os.path.join(log_dir, f"spatialsnake_{self.runid}_{self._log_label()}_{timestamp}_runlog.log")
        return self.logname

    def append_to_log(self, text: str):
        """Append text to the current run log without changing console behavior."""
        if not self.log_enabled:
            return
        try:
            logname = self._ensure_logname()
            os.makedirs(os.path.dirname(logname), exist_ok=True)
            with open(logname, "a", encoding="utf-8") as f:
                f.write(text)
        except Exception as e:
            logger.error(f"Failed to write log file: {e}")

    def append_log_rule(self, char: str = "="):
        self.append_to_log(char * LOG_WIDTH + "\n")

    def append_log_kv(self, key: str, value: Any):
        self.append_to_log(f"  {key.ljust(18)}: {value}\n")

    def format_datetime(self, value):
        if value is None:
            return "NA"
        return value.strftime("%Y-%m-%d %H:%M:%S")

    def format_peak_memory(self):
        if self.resource_after is None:
            return "NA"
        peak = float(getattr(self.resource_after, "ru_maxrss", 0))
        if peak <= 0:
            return "NA"
        # Linux reports ru_maxrss in KiB; macOS reports bytes.
        peak_mb = peak / (1024 * 1024) if sys.platform == "darwin" else peak / 1024
        if peak_mb >= 1024:
            return f"{peak_mb / 1024:.2f} GB"
        return f"{peak_mb:.2f} MB"

    def format_cpu_time(self, attr: str):
        if self.resource_before is None or self.resource_after is None:
            return "NA"
        before = float(getattr(self.resource_before, attr, 0))
        after = float(getattr(self.resource_after, attr, 0))
        return f"{max(0.0, after - before):.2f} sec"

    def write_log_header(self):
        """Write run metadata before subprocess output starts."""
        if not self.log_enabled or self.log_header_written:
            return
        if self.started_at is None:
            self.started_at = datetime.datetime.now()

        option_val = self.arguments.get("--option", "unknown")
        channel_val = self.arguments.get("<command>") or ("useful_tool" if self.arguments.get("useful_tool") else "unknown")
        run_type_val = self.arguments.get("<TYPE>", "unknown")
        results_folder = self.parameters.get("results_folder", self.arguments.get("--results_folder", "results"))
        python_version = sys.version.replace("\n", " ")
        configfile = self.configfile_path or self.arguments.get("--configfile", "")

        self.append_log_rule("=")
        self.append_to_log("SPATIALSNAKE RUN LOG\n")
        self.append_log_rule("=")
        self.append_to_log("\nRun information\n")
        self.append_log_kv("Run ID", self.runid)
        self.append_log_kv("Version", __version__)
        self.append_log_kv("Started at", self.format_datetime(self.started_at))
        self.append_log_kv("Working directory", os.getcwd())
        self.append_log_kv("Python version", python_version)
        self.append_to_log("\nAnalysis context\n")
        self.append_log_kv("Command type", channel_val)
        self.append_log_kv("Analysis option", option_val)
        self.append_log_kv("Data type", run_type_val)
        self.append_log_kv("Config file", configfile)
        self.append_log_kv("Results folder", results_folder)
        self.append_to_log("\nCommand\n")
        self.append_log_kv("CLI command", " ".join(sys.argv))
        self.append_log_kv("Backend command", self.cmd_str)
        self.append_to_log("\nParameters\n")
        for key, value in sorted(self.parameters.items()):
            if value is not None and value != "":
                self.append_log_kv(str(key), value)
        self.append_to_log("\n")
        self.append_log_rule("-")
        self.append_to_log("Pipeline output\n")
        self.append_log_rule("-")
        self.log_header_written = True

    def write_to_log(self, start_time: float):
        """Write execution summary to log file."""
        if not self.log_enabled:
            return

        if not self.log_header_written:
            self.write_log_header()

        stop_time = timeit.default_timer()
        elapsed_mins = (stop_time - start_time) / 60
        if self.finished_at is None:
            self.finished_at = datetime.datetime.now()

        self.append_to_log("\n")
        self.append_log_rule("-")
        self.append_to_log("Run summary\n")
        self.append_log_rule("-")
        self.append_log_kv("Status", self.execution_status)
        self.append_log_kv("Return code", self.return_code if self.return_code is not None else "NA")
        if self.error_message:
            self.append_log_kv("Error message", self.error_message)
        self.append_log_kv("Finished at", self.format_datetime(self.finished_at))
        self.append_log_kv("Elapsed time", f"{elapsed_mins:.2f} mins")
        self.append_log_kv("Peak memory", self.format_peak_memory())
        self.append_log_kv("CPU user time", self.format_cpu_time("ru_utime"))
        self.append_log_kv("CPU system time", self.format_cpu_time("ru_stime"))
        if self.execution_status == "FAILED":
            self.append_to_log("\nFailure hint\n")
            self.append_to_log("  The command failed after the pipeline output shown above.\n")
            self.append_to_log("  Please check the last ERROR, Traceback, or Snakemake rule message in this log.\n")
            self.append_log_kv("Failed command", self.cmd_str)
        self.append_to_log("\nLog file\n")
        self.append_log_kv("Path", self._ensure_logname())
        self.append_log_rule("=")

    def run_subprocess(self):
        """Run prepared command while teeing stdout/stderr into the run log."""
        process = subprocess.Popen(
            str(self.cmd_str),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )

        if process.stdout:
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                self.append_to_log(line)

        return process.wait()

    def execute(self):
        """Execute the prepared command."""
        start_time = timeit.default_timer()
        exit_code = 0
        try:
            self.started_at = datetime.datetime.now()
            if resource is not None:
                self.resource_before = resource.getrusage(resource.RUSAGE_CHILDREN)
            self.prepare_arguments()
            self.write_log_header()
            logger.info(f"Executing command: {self.cmd_str}")
            self.append_to_log(f"Executing command: {self.cmd_str}\n")
            exit_code = self.run_subprocess()
            self.return_code = exit_code
            if exit_code == 0:
                self.execution_status = "SUCCESS"
                logger.info("Command execution finished successfully")
            else:
                self.execution_status = "FAILED"
                self.error_message = f"Command execution failed with return code {exit_code}"
                logger.error(self.error_message)
        except Exception as e:
            exit_code = 1
            self.return_code = exit_code
            self.execution_status = "FAILED"
            self.error_message = str(e)
            if isinstance(e, FileNotFoundError):
                logger.error(str(e))
            else:
                logger.error(f"An unexpected error occurred: {e}")
        finally:
            self.finished_at = datetime.datetime.now()
            if resource is not None:
                self.resource_after = resource.getrusage(resource.RUSAGE_CHILDREN)
            self.write_to_log(start_time)
        if exit_code != 0:
            sys.exit(exit_code)


class WorkflowRunner(BaseCommandLine):
    """Handles Snakemake workflow execution."""

    def __init__(self, arguments: Dict[str, Any]):
        super().__init__(arguments)
        self.snakemake_cmd = "snakemake --rerun-incomplete -k "
        self.cmd_str = self.snakemake_cmd

    def add_config_argument(self):
        self.cmd_str += " --config " + " ".join(shlex.quote(str(item)) for item in self.config)

    def load_configfile(self):
        step = self.arguments.get("--option")
        configfile = None
        requested_configfile = self.arguments.get("--configfile")
        if cli_option_provided("--configfile"):
            if not requested_configfile or not os.path.isfile(requested_configfile):
                raise FileNotFoundError(f"Config file not found: {requested_configfile}")
            configfile = requested_configfile
            self.cmd_str += f" --configfile={configfile}"
            self.configfile_loaded = True
        else:
            # Fallback to default env config
            default_config = os.path.join(SPATIALSNAKE_PATH, f"workflow/envs/{step}.yaml")
            self.cmd_str += f" --configfile={default_config}"
            configfile = default_config
            self.arguments["--configfile"] = default_config
        self.configfile_path = configfile
            
        if configfile and os.path.exists(configfile):
            try:
                with open(configfile) as f:
                    self.parameters = yaml.load(f, Loader=SafeLoader) or {}
            except Exception as e:
                logger.warning(f"Failed to load config file {configfile}: {e}")

    def prepare_arguments(self):
        jobs = self.arguments.get('--jobs') or 16
        self.cmd_str += f" -j {jobs} "
        self.cmd_str += f" -s {os.path.join(SPATIALSNAKE_PATH, 'workflow/Snakefile')} "
        
        self.load_configfile()
        
        if self.arguments.get('--option') in ["integrate", "preprocess", "clustering", "reclustering", "annotation_help", "annotation", "compare_stage", "all", "advance_analysis"]:
            self.config.append(f"sample_list={self.arguments.get('<INPUT_FILE>')}")
            self.parameters["sample_list"] = self.arguments.get("<INPUT_FILE>")
            
        spatialsnake_path = f"{SPATIALSNAKE_PATH}/"
        self.config.append(f"spatialsnake_path={spatialsnake_path}")
        self.parameters["spatialsnake_path"] = spatialsnake_path
        
        # Parse dynamic arguments. Only user-provided CLI options override yaml values.
        exclude_keys = {
            "--jobs", "--configfile", "--option", "--unlock", "--remove", 
            "--dry", "--dry-run", "--extended", "--torch-backend", "--help", "--version", "<INPUT_FILE>", "<command>",
            "--install-packages", "<TYPE>", "useful_tool", "<INTEGRATED_FILE>"
        }
        for key, value in self.arguments.items():
            if key in exclude_keys:
                continue
                
            clean_key = config_key_from_cli(key)
            user_configfile = self.configfile_loaded and cli_option_provided("--configfile")
            explicit_cli_value = value is not None and cli_option_provided(key)

            # Default templates keep module-specific filtering. A user configfile
            # may be partial, so an explicitly provided CLI value can add a key.
            if clean_key not in self.parameters and not (user_configfile and explicit_cli_value):
                continue

            if not explicit_cli_value:
                continue

            self.config.append(f"{clean_key}={value}")
            self.parameters[clean_key] = str(value)
        self.config.append(f"runid={self.runid}")
        self.parameters["runid"] = self.runid
        
        if self.arguments.get("--option"):
            self.config.append(f"option={self.arguments['--option']}")
            self.parameters["option"] = self.arguments["--option"]
            
        type_value = self.arguments.get('<TYPE>')
        type_alias = {
            "merfish": "Merfish",
            "merscope": "Merfish",
            "MERFISH": "Merfish",
            "MERSCOPE": "Merfish",
            "StereoSeq": "stereoseq",
            "Stereo-seq": "stereoseq",
            "stereo-seq": "stereoseq",
        }
        type_norm = type_alias.get(type_value, type_value)
        self.config.append(f"channel={self.arguments.get('<command>')}")
        self.config.append(f"run_type={type_norm}")
        self.parameters["channel"] = self.arguments.get("<command>")
        self.parameters["run_type"] = type_norm

        validation_config = dict(self.parameters)
        validation_config["jobs"] = jobs
        try:
            validate_workflow_config(validation_config)
        except ConfigValidationError as exc:
            logger.error(str(exc))
            sys.exit(1)
        
        if self.arguments.get("--dry"):
            self.cmd_str += " -n "
            self.log_enabled = False
        if self.arguments.get("--unlock"):
            self.cmd_str += " --unlock "
            self.log_enabled = False
        if self.arguments.get("--remove"):
            self.cmd_str += " --delete-all-output "
            self.log_enabled = False
            
        self.add_config_argument()


class ToolRunner(BaseCommandLine):
    """Handles useful_tool execution."""

    OPTIONAL_EMPTY_TOOL_PARAMS = {
        "annotation_cols",
        "annotation_csv",
        "barcodes",
        "coordinate_system",
        "neighbor_col",
        "neighbor_types",
        "roi_csv",
        "roi_label_col",
        "roi_region",
        "sample_ids",
        "label_cols",
        "table_key",
    }

    def __init__(self, arguments: Dict[str, Any]):
        super().__init__(arguments)
        self.cmd_str = "python "

    def add_config_argument(self):
        self.cmd_str += " ".join(self.config)

    def load_configfile(self):
        tool = self.arguments.get("--option")
        if tool:
            self.config.append(os.path.join(SPATIALSNAKE_PATH, f"workflow/function/{tool}.py"))
            
        configfile = None
        requested_configfile = self.arguments.get("--configfile")
        if cli_option_provided("--configfile"):
            if not requested_configfile or not os.path.isfile(requested_configfile):
                raise FileNotFoundError(f"Config file not found: {requested_configfile}")
            configfile = requested_configfile
            self.configfile_loaded = True
            logger.info(f"Using config file: {configfile}")
        else:
            default_config = os.path.join(SPATIALSNAKE_PATH, f"workflow/envs/{tool}.yaml")
            configfile = default_config
            self.arguments["--configfile"] = default_config
            logger.info(f"Using default config file: {configfile}")
        self.configfile_path = configfile
            
        if configfile and os.path.exists(configfile):
            try:
                with open(configfile) as f:
                    self.parameters = yaml.load(f, Loader=SafeLoader) or {}
            except Exception as e:
                logger.warning(f"Failed to load config file {configfile}: {e}")

    def build_subprocess_cmd(self):
        if self.arguments.get("--option") in ["merge", "transform"]:
            self.config.extend(['--INPUT'])
            inputs = self.arguments.get('<INPUT>', [])
            if isinstance(inputs, list):
                self.config.extend(shlex.quote(str(value)) for value in inputs)
            else:
                 self.config.append(shlex.quote(str(inputs)))
        elif self.arguments.get("--option") == "splitting":
            inputs = self.arguments.get('<INPUT>')
            input_value = inputs[0] if isinstance(inputs, list) and len(inputs) > 0 else inputs
            self.config.append(f"--INPUT_FILE {shlex.quote(str(input_value))}")
        else:
            # Original code: cmd.append("--INPUT {}".format(arguments['<INPUT>']))
            # <INPUT> is usually a list in docopt if ... is used, but here for useful_tool it is <INPUT>...
            # But in the original code for non-merge/transform: 
            # cmd.append("--INPUT {}".format(arguments['<INPUT>']))
            # If <INPUT> is a list, format might behave weirdly if not handled.
            # Assuming <INPUT> is a list, let's join it or take first? 
            # The original code used arguments['<INPUT>'] directly in format.
            # If docopt returns list for <INPUT>..., then format will stringify the list.
            inputs = self.arguments.get('<INPUT>')
            self.config.append(f"--INPUT {inputs}")

    def prepare_arguments(self):
        self.load_configfile()
        self.build_subprocess_cmd()
        
        exclude_keys = {"--jobs", "--configfile", "--option", "--dry-run", "--extended", "--torch-backend", "useful_tool", "<INPUT>"}

        cli_config_keys = []
        for key, value in self.arguments.items():
            if key in exclude_keys:
                continue

            clean_key = config_key_from_cli(key)
            user_configfile = self.configfile_loaded and cli_option_provided("--configfile")
            explicit_cli_value = value is not None and cli_option_provided(key)
            if clean_key not in self.parameters and not (user_configfile and explicit_cli_value):
                continue

            if clean_key not in cli_config_keys:
                cli_config_keys.append(clean_key)

            if explicit_cli_value:
                self.parameters[clean_key] = str(value)

        self.validate_effective_tool_config()

        for clean_key in cli_config_keys:
            value = self.parameters.get(clean_key)
            if self.should_skip_tool_param(clean_key, value):
                continue
            self.config.append(f"--{clean_key} {self.format_tool_value(value)}")

        self.add_config_argument()

    def should_skip_tool_param(self, key: str, value: Any) -> bool:
        if key in {"INPUT", "INPUT_FILE"}:
            return True
        if value is None:
            return True
        if key in self.OPTIONAL_EMPTY_TOOL_PARAMS and str(value).strip().lower() in {"", "none", "null"}:
            return True
        return False

    def validate_effective_tool_config(self):
        option = self.arguments.get("--option")
        if option == "splitting":
            try:
                validate_splitting_config(self.parameters)
            except ConfigValidationError as exc:
                logger.error(str(exc))
                sys.exit(1)
        if option == "merge":
            try:
                validate_merge_config(self.parameters)
            except ConfigValidationError as exc:
                logger.error(str(exc))
                sys.exit(1)
        if option == "transform":
            try:
                validate_transform_config(self.parameters)
            except ConfigValidationError as exc:
                logger.error(str(exc))
                sys.exit(1)

    def format_tool_value(self, value: Any) -> str:
        return shlex.quote(str(value))


def validate_workflow_arguments(arguments: Dict[str, Any]) -> bool:
    """Validate arguments for workflow commands."""
    input_file = arguments.get("<INPUT_FILE>")
    if not input_file or not os.path.isfile(input_file):
        logger.error("未提供必要样本信息，请提供 sample.txt 或通过命令指定样本表路径")
        if input_file:
            logger.error(f"Sample list file not found: {input_file}")
        return False

    option = arguments.get("--option")
    type_arg = arguments.get("<TYPE>")
    normalized_type = str(type_arg).strip().lower().replace("-", "_")
    unsupported_type = "slide" + "_seq"
    if normalized_type in {unsupported_type, "slide" + "seq"}:
        logger.error("Slide" + "-seq input is not supported by this version of Spatialsnake.")
        return False
    type_alias = {
        "merfish": "Merfish",
        "merscope": "Merfish",
        "MERFISH": "Merfish",
        "MERSCOPE": "Merfish",
        "StereoSeq": "stereoseq",
        "Stereo-seq": "stereoseq",
        "stereo-seq": "stereoseq",
    }
    type_arg = type_alias.get(type_arg, type_arg)
    arguments["<TYPE>"] = type_arg
    valid_types = ['visium', 'visium_segment', 'visium_HD', 'xenium', 'Merfish', 'stereoseq']
    if type_arg not in valid_types:
        logger.error(f"Invalid spatialdata type. Valid types: {', '.join(valid_types)}")
        return False

    if option and option not in VALID_OPTIONS:
        logger.error("Invalid option selected.")
        logger.info(f"Correct options include: {' '.join(VALID_OPTIONS)}")
        return False

    return True


def validate_tool_arguments(arguments: Dict[str, Any]) -> bool:
    """Validate arguments for useful_tool."""
    # Note: <INPUT> is a list for 'useful_tool' command in docopt because of <INPUT>...
    # But check_arguments_inputfile in original code treated it as single path in one check: os.path.exists(arguments["<INPUT>"])
    # If <INPUT> is a list, os.path.exists will fail. 
    # Let's check how docopt parses `spatialsnake useful_tool ... <INPUT>...`
    # It returns a list.
    # The original code:
    # if not os.path.exists(arguments["<INPUT>"]):
    # This implies arguments["<INPUT>"] was expected to be a string or the original code was buggy for multiple inputs?
    # Or maybe <INPUT>... means list, but if user provides one, it is a list of one.
    # I will iterate if it is a list.
    
    inputs = arguments.get("<INPUT>")
    if isinstance(inputs, list):
        for inp in inputs:
            if not os.path.exists(inp):
                logger.error(f"Input file not found: {inp}")
                return False
    elif isinstance(inputs, str):
        if not os.path.exists(inputs):
            logger.error(f"Input file not found: {inputs}")
            return False

    option = arguments.get("--option")
    valid_tool_options = ["splitting", "transform", "merge"]
    if option not in valid_tool_options:
        logger.error(f"Invalid option for useful_tool. Valid options: {', '.join(valid_tool_options)}")
        return False

    return True


INSTALL_PACKAGES_CONFIG = os.path.join(SPATIALSNAKE_PATH, "workflow/envs/install_packages.yaml")


def conda_package_satisfied(name: str, version: Optional[str]) -> bool:
    meta_dir = os.path.join(sys.prefix, "conda-meta")
    if name == "gfortran_linux-64":
        compiler = os.path.join(sys.prefix, "bin", "x86_64-conda-linux-gnu-gfortran")
        if not os.path.isfile(compiler):
            return False
        for filename in os.listdir(meta_dir):
            if filename.startswith("gfortran_linux-64-") and filename.endswith(".json"):
                try:
                    with open(os.path.join(meta_dir, filename), encoding="utf-8") as f:
                        return json.load(f).get("version") == version
                except (OSError, ValueError):
                    return False
        return False
    if name.startswith(("r-", "bioconductor-")):
        if not os.path.isdir(meta_dir):
            return False
        for filename in os.listdir(meta_dir):
            if not filename.endswith(".json"):
                continue
            try:
                with open(os.path.join(meta_dir, filename), encoding="utf-8") as f:
                    package = json.load(f)
            except (OSError, ValueError):
                continue
            if package.get("name") == name:
                return version is None or package.get("version") == version
        return False
    try:
        installed_version = importlib.metadata.version(name)
        return version is None or installed_version == version
    except importlib.metadata.PackageNotFoundError:
        return False


def install_conda_requirements(packages: Dict[str, Optional[str]], label: str, dry_run: bool = False):
    """Repair compatibility packages that must remain Conda-managed."""
    pending = [f"{name}={version}" if version else name for name, version in packages.items()
               if not conda_package_satisfied(name, version)]
    if not pending:
        logger.info(f"All {label} are already installed.")
        return

    conda = os.environ.get("CONDA_EXE")
    if not conda or not os.path.isfile(conda):
        conda = os.path.join(
            os.path.dirname(os.path.dirname(sys.prefix)), "bin", "conda"
        )
    if not os.path.isfile(conda):
        conda = shutil.which("conda")
    if not conda or not os.path.isfile(conda):
        logger.error("Conda executable not found; activate the Bioconda environment and retry.")
        sys.exit(1)
    cmd = [conda, "install", "-p", sys.prefix, "--freeze-installed",
           "--use-index-cache",
           "--override-channels", "-c", "conda-forge", "-c", "bioconda",
           *pending, "-y"]
    if dry_run:
        logger.info(f"DRY RUN: would install {label}: {' '.join(cmd)}")
        return
    logger.info(f"Installing {label}: {', '.join(pending)}")
    try:
        subprocess.check_call(cmd, env=get_subprocess_env())
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install {label}: {e}")
        sys.exit(e.returncode)


def check_conda_runtime(dry_run: bool = False):
    if sys.version_info[:2] != (3, 12):
        logger.error("This install profile is validated for Python 3.12. "
                     "Create a fresh Bioconda environment with python=3.12.11.")
        sys.exit(2)
    if dry_run:
        return
    try:
        subprocess.check_call(
            [get_rscript_executable(), "-e",
             'stopifnot(getRversion() >= "4.4.0", getRversion() < "4.5.0")'],
            env=get_subprocess_env(), stdout=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        logger.error("This install profile needs R 4.4; create a fresh Conda "
                     "environment with r-base=4.4.0.")
        sys.exit(2)


def repair_numpy_runtime(constraints: List[str], dry_run: bool = False):
    target = next((item for item in constraints if item.startswith("numpy==")), None)
    if not target:
        return
    if dry_run:
        logger.info(f"DRY RUN: would verify runtime {target} after Conda/pip installation")
        return
    try:
        actual = subprocess.check_output(
            [sys.executable, "-c", "import numpy; print(numpy.__version__)"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        actual = None
    if actual == target.split("==", 1)[1]:
        return
    logger.warning(f"NumPy runtime is {actual or 'unavailable'}; repairing to {target}")
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "--force-reinstall", "--no-deps", target]
    )


def get_torch_runtime() -> Tuple[str, bool]:
    """Return the installed Torch CUDA runtime and whether CUDA is usable."""
    try:
        output = subprocess.check_output(
            [
                sys.executable, "-c",
                "import torch; print(torch.version.cuda or ''); "
                "print('1' if torch.cuda.is_available() else '0')",
            ],
            text=True,
            stderr=subprocess.DEVNULL,
        ).splitlines()
    except (OSError, subprocess.CalledProcessError):
        return "", False
    runtime = output[0].strip() if output else ""
    available = len(output) > 1 and output[1].strip() == "1"
    return runtime, available


def nvidia_gpu_available() -> bool:
    """Check for an accessible NVIDIA GPU without requiring PyTorch."""
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return False
    try:
        output = subprocess.check_output(
            [nvidia_smi, "-L"], text=True, stderr=subprocess.DEVNULL, timeout=10,
        )
        return bool(output.strip())
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


def resolve_torch_backend(requested: str) -> str:
    backend = str(requested or "auto").strip().lower()
    if backend not in {"auto", "cpu", "cu128"}:
        logger.error("--torch-backend must be one of: auto, cpu, cu128")
        sys.exit(2)
    if backend != "auto":
        return backend

    _, torch_cuda_available = get_torch_runtime()
    selected = "cu128" if torch_cuda_available or nvidia_gpu_available() else "cpu"
    logger.info(f"Auto-selected PyTorch backend: {selected}")
    return selected


def install_pytorch(
    config: Dict[str, Any], backend: str = "auto", dry_run: bool = False
) -> str:
    """Install the tested CPU or CUDA PyTorch build before dependent packages."""
    requirement = config.get("requirement", "torch==2.8.0")
    if backend == "auto":
        if package_satisfies_requirement(requirement):
            cuda_runtime, _ = get_torch_runtime()
            logger.info(f"Reusing user-installed {requirement}; CUDA runtime: {cuda_runtime or 'CPU'}")
            return "cu128" if cuda_runtime else "cpu"
        logger.error(
            f"Install {requirement} first, choosing --index-url "
            "https://download.pytorch.org/whl/cpu or "
            "https://download.pytorch.org/whl/cu128."
        )
        if dry_run:
            return "manual"
        sys.exit(2)
    selected = resolve_torch_backend(backend)
    indexes = config.get("index_urls", {})
    index_url = indexes.get(
        selected, f"https://download.pytorch.org/whl/{selected}"
    )
    cuda_runtime, _ = get_torch_runtime()
    if package_satisfies_requirement(requirement):
        matching_build = (selected == "cpu" and not cuda_runtime) or (
            selected == "cu128" and cuda_runtime.startswith("12.8")
        )
        if matching_build:
            build = "CPU" if selected == "cpu" else f"CUDA {cuda_runtime}"
            logger.info(f"{requirement} is already installed with the requested {build} build.")
            return selected

    cmd = [
        sys.executable, "-m", "pip", "install", "--force-reinstall",
        requirement, "--index-url", index_url,
        "--upgrade-strategy", "only-if-needed",
    ]
    if dry_run:
        logger.info(f"DRY RUN: would install {selected} PyTorch: {' '.join(cmd)}")
        return selected
    logger.info(f"Installing {selected} PyTorch from {index_url}: {requirement}")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install {selected} PyTorch: {e}")
        sys.exit(e.returncode)

    cuda_runtime, _ = get_torch_runtime()
    matching_build = (selected == "cpu" and not cuda_runtime) or (
        selected == "cu128" and cuda_runtime.startswith("12.8")
    )
    if not package_satisfies_requirement(requirement) or not matching_build:
        logger.error(f"Installed PyTorch does not match the requested {selected} build.")
        sys.exit(1)
    return selected


def install_packages(
    extended: bool = False, torch_backend: str = "auto", dry_run: bool = False
):
    """Complete dependencies using the installation route for this package."""
    torch_backend = str(torch_backend or "auto").strip().lower()
    if torch_backend not in {"auto", "cpu", "cu128"}:
        logger.error("--torch-backend must be one of: auto, cpu, cu128")
        sys.exit(2)
    origin = detect_spatialsnake_install_origin()
    logger.info(f"Detected Spatialsnake install origin: {origin}")
    config = load_install_packages_config()
    conda_config = config.get("conda", {})
    check_conda_runtime(dry_run=dry_run)

    if origin == "conda":
        pip_config = config.get("pip", {})
        constraints = config.get("constraints", [])

        conda_packages = dict(conda_config.get("compatibility", {}))
        if extended:
            conda_packages.update(conda_config.get("extended", {}))
        install_conda_requirements(
            conda_packages, label="Conda compatibility and R packages",
            dry_run=dry_run,
        )
        core_requirements = list(pip_config.get("conda_minimal", []))
        for constraint in constraints:
            name, _ = parse_pinned_requirement(constraint)
            if name == "torch":
                continue
            try:
                importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                continue
            if not package_satisfies_requirement(constraint) and constraint not in core_requirements:
                core_requirements.append(constraint)
        install_pip_requirements(
            core_requirements,
            label="core Python packages", dry_run=dry_run, constraints=constraints,
        )
        install_pip_requirements(
            pip_config.get("conda_force_pypi", []),
            label="validated Numba/llvmlite wheels", dry_run=dry_run,
            constraints=constraints, force_reinstall=True, no_deps=True,
            require_pip_installer=True,
        )
        repair_numpy_runtime(constraints, dry_run=dry_run)
        run_conda_r_core_install_script(dry_run=dry_run)
        if not extended:
            logger.info("Core analysis packages handled. Use --extended for "
                        "cell2location, downstream Python tools, RCTD and CellChat.")
            return

        install_pytorch(
            pip_config.get("pytorch", {}), backend=torch_backend, dry_run=dry_run
        )
        install_pip_requirements(
            pip_config.get("conda_extended", []),
            label="extended Python packages", dry_run=dry_run,
            constraints=constraints,
        )
        repair_numpy_runtime(constraints, dry_run=dry_run)
        _install_pybanksy_if_needed(
            requirement=config.get("pybanksy", {}).get("requirement", "pybanksy==1.3.5"),
            dry_run=dry_run,
        )
        run_conda_r_install_script(dry_run=dry_run)
        run_conda_r_github_install_script(dry_run=dry_run)
        return

    if extended:
        logger.info(
            "--extended does not install Python extras for PyPI/source environments. "
            "Use `python -m pip install \"spatialsnake[extended]\"` before this command."
        )
    _install_pybanksy_if_needed(
        requirement=config.get("pybanksy", {}).get("requirement", "pybanksy==1.3.5"),
        dry_run=dry_run,
    )
    run_r_install_script(dry_run=dry_run)


def detect_spatialsnake_install_origin() -> str:
    meta_dir = os.path.join(sys.prefix, "conda-meta")
    if not os.path.isdir(meta_dir):
        return "pypi"

    for filename in os.listdir(meta_dir):
        if not filename.startswith("spatialsnake-") or not filename.endswith(".json"):
            continue
        path = os.path.join(meta_dir, filename)
        try:
            with open(path, "r", encoding="utf-8") as f:
                metadata = json.load(f)
        except Exception:
            logger.warning(f"Could not read conda metadata file: {path}")
            continue

        source = " ".join(str(metadata.get(key, "")) for key in ("channel", "platform", "base_url", "url")).lower()
        if "pypi" not in source:
            return "conda"
    return "pypi"


def load_install_packages_config() -> Dict[str, Any]:
    try:
        with open(INSTALL_PACKAGES_CONFIG, "r", encoding="utf-8") as f:
            return yaml.load(f, Loader=SafeLoader) or {}
    except FileNotFoundError:
        logger.error(f"Install package config not found: {INSTALL_PACKAGES_CONFIG}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to load install package config: {e}")
        sys.exit(1)


def parse_pinned_requirement(requirement: str):
    if "==" not in requirement:
        return requirement, None
    package_name, version = requirement.split("==", 1)
    return package_name.strip(), version.strip()


def package_satisfies_requirement(requirement: str) -> bool:
    package_name, version = parse_pinned_requirement(requirement)
    if version is None:
        return importlib.util.find_spec(package_name) is not None
    try:
        installed_version = importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return False
    return installed_version.split("+", 1)[0] == version


def package_installer(package_name: str) -> str:
    try:
        distribution = importlib.metadata.distribution(package_name)
    except importlib.metadata.PackageNotFoundError:
        return ""
    try:
        return (distribution.read_text("INSTALLER") or "").strip().lower()
    except Exception:
        return ""


def package_satisfies_pip_requirement(requirement: str) -> bool:
    package_name, _ = parse_pinned_requirement(requirement)
    return package_satisfies_requirement(requirement) and package_installer(package_name) == "pip"


def install_pip_requirements(
    requirements: List[str],
    label: str,
    dry_run: bool = False,
    constraints: List[str] = None,
    force_reinstall: bool = False,
    no_deps: bool = False,
    require_pip_installer: bool = False,
    extra_index_url: Optional[str] = None,
):
    satisfies = package_satisfies_pip_requirement if require_pip_installer else package_satisfies_requirement
    pending = [req for req in requirements if not satisfies(req)]
    if not pending:
        logger.info(f"All {label} are already installed with the required versions.")
        return

    cmd = [sys.executable, "-m", "pip", "install"]
    if force_reinstall:
        cmd.append("--force-reinstall")
    if no_deps:
        cmd.append("--no-deps")
    if extra_index_url:
        cmd.extend(["--extra-index-url", extra_index_url])
    cmd.extend(pending)
    cmd.extend(["--upgrade-strategy", "only-if-needed"])
    if dry_run:
        logger.info(f"DRY RUN: would install {label}: {' '.join(cmd)}")
        if constraints:
            logger.info(f"DRY RUN: would constrain pip resolver with {len(constraints)} pinned versions from install_packages.yaml")
        return

    constraint_path = None
    if constraints:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False) as f:
            f.write("\n".join(constraints))
            f.write("\n")
            constraint_path = f.name
        cmd.extend(["--constraint", constraint_path])

    logger.info(f"Installing {label}: {', '.join(pending)}")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install {label}: {e}")
        sys.exit(e.returncode)
    finally:
        if constraint_path:
            try:
                os.unlink(constraint_path)
            except OSError:
                pass


def run_r_install_script(dry_run: bool = False):
    r_script_path = os.path.join(SPATIALSNAKE_PATH, "workflow/scripts/install_packages.R")
    rscript = get_rscript_executable()
    env = get_subprocess_env()
    if dry_run:
        logger.info(f"DRY RUN: would run R script with {rscript}: {r_script_path}")
        return
    try:
        logger.info(f"Running R script with {rscript}: {r_script_path}")
        subprocess.check_call([rscript, r_script_path], env=env)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install R packages: {e}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Rscript not found or script missing.")
        sys.exit(1)


def run_conda_r_core_install_script(dry_run: bool = False):
    r_script_path = os.path.join(
        SPATIALSNAKE_PATH, "workflow/scripts/install_packages_core_conda.R"
    )
    rscript = get_rscript_executable()
    if dry_run:
        logger.info(f"DRY RUN: would run core R script with {rscript}: {r_script_path}")
        return
    try:
        logger.info(f"Running core R script with {rscript}: {r_script_path}")
        subprocess.check_call([rscript, r_script_path], env=get_subprocess_env())
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install core R packages: {e}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Rscript not found or core R script missing.")
        sys.exit(1)


def run_conda_r_install_script(dry_run: bool = False):
    r_script_path = os.path.join(SPATIALSNAKE_PATH, "workflow/scripts/install_packages_conda.R")
    rscript = get_rscript_executable()
    env = get_subprocess_env()
    if dry_run:
        logger.info(f"DRY RUN: would run conda-only R script with {rscript}: {r_script_path}")
        return
    try:
        logger.info(f"Running conda-only R script with {rscript}: {r_script_path}")
        subprocess.check_call([rscript, r_script_path], env=env)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install conda-only R packages: {e}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Rscript not found or conda-only R script missing.")
        sys.exit(1)


def run_conda_r_github_install_script(dry_run: bool = False):
    r_script_path = os.path.join(SPATIALSNAKE_PATH, "workflow/scripts/install_packages_github_conda.R")
    rscript = get_rscript_executable()
    env = get_subprocess_env()
    if dry_run:
        logger.info(f"DRY RUN: would run conda-only R/GitHub script with {rscript}: {r_script_path}")
        return
    try:
        logger.info(f"Running conda-only R/GitHub script with {rscript}: {r_script_path}")
        subprocess.check_call([rscript, r_script_path], env=env)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install conda-only R/GitHub packages: {e}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Rscript not found or conda-only R/GitHub script missing.")
        sys.exit(1)


def get_rscript_executable() -> str:
    env_rscript = os.path.join(sys.prefix, "bin", "Rscript")
    if os.path.exists(env_rscript):
        return env_rscript
    return "Rscript"


def get_subprocess_env() -> Dict[str, str]:
    env = os.environ.copy()
    env_bin = os.path.join(sys.prefix, "bin")
    if os.path.isdir(env_bin):
        env["PATH"] = env_bin + os.pathsep + env.get("PATH", "")
        env.setdefault("CONDA_PREFIX", sys.prefix)
    return env


def _install_pybanksy_if_needed(
    requirement: str = "pybanksy==1.3.5", dry_run: bool = False
):
    """Install pybanksy without resolver changes; its dependencies are already
    satisfied by the extended environment.

    The distribution is named pybanksy, but its import package is `banksy`.
    Check distribution metadata so repeated runs are idempotent.
    """
    _, expected_version = parse_pinned_requirement(requirement)
    try:
        if expected_version and importlib.metadata.version("pybanksy") == expected_version:
            return
    except importlib.metadata.PackageNotFoundError:
        pass

    if dry_run:
        logger.info(f"DRY RUN: would install pybanksy with: {sys.executable} -m pip install {requirement} --no-deps")
        return

    logger.info("Installing pybanksy (BANKSY spatial clustering) ...")
    try:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", requirement, "--no-deps"],
        )
        logger.info("pybanksy installed successfully.")
    except subprocess.CalledProcessError:
        logger.warning(
            "Failed to install pybanksy automatically. "
            f"You can install it manually later: pip install {requirement} --no-deps"
        )


def generate_config_file(arguments: Dict[str, Any]):
    """Generate configuration file."""
    step = arguments.get("--option")
    command = arguments.get("<command>")
    
    # Logic from original main() to determine 'step'
    if arguments.get("produce-file"):
         if command not in ['useful_tool', 'transform']:
             # If command is not useful_tool/transform, use --option as step.
             # But arguments['<command>'] is None if running `spatialsnake produce-file ...`
             # The usage for produce-file is: `spatialsnake produce-file [--option=<analysis_option>]`
             pass

    # Simplified logic:
    if step not in VALID_OPTIONS + ["all"]:
        logger.error("Please set correct params: --option=<step_name> or --option=all")
        return False

    logger.info(f"Generating config.yaml file for: {step}...")
    
    try:
        if step == "all":
            src = os.path.join(SPATIALSNAKE_PATH, "config.yaml")
            dst = "config.yaml"
            shutil.copyfile(src, dst)
        else:
            src = os.path.join(SPATIALSNAKE_PATH, f"workflow/envs/{step}.yaml")
            dst = f"{step}.yaml"
            shutil.copyfile(src, dst)
            
        logger.info("You can use this as a config-file for a spatialsnake run.")
        logger.info("***How to set your own params:***")
        logger.info("Add params: --configfile <file-path> in the command line")
        return True
        
    except FileNotFoundError:
        logger.error(f"Source config file not found: {src}")
        return False
    except Exception as e:
        logger.error(f"Error generating config file: {e}")
        return False


def main():
    """Main entry point."""
    try:
        cli_arguments = docopt(
            __doc__,
            argv=normalize_cli_equals_tokens(sys.argv[1:]),
            version=__version__,
        )
    except Exception as e:
        logger.error(f"Error parsing arguments: {e}")
        sys.exit(2)

    cli_arguments["--option"] = OPTION_ALIASES.get(cli_arguments.get("--option"), cli_arguments.get("--option"))
    cli_arguments["--runpipe"] = RUNPIPE_ALIASES.get(cli_arguments.get("--runpipe"), cli_arguments.get("--runpipe"))

    if cli_arguments.get("produce-file"):
        if not generate_config_file(cli_arguments):
            sys.exit(2)
        return

    if cli_arguments.get("install-packages"):
        install_packages(
            extended=bool(cli_arguments.get("--extended")),
            torch_backend=cli_arguments.get("--torch-backend") or "auto",
            dry_run=bool(cli_arguments.get("--dry-run")),
        )
        return

    if cli_arguments.get("useful_tool"):
        if not validate_tool_arguments(cli_arguments):
             sys.exit(2)
        runner = ToolRunner(cli_arguments)
        runner.execute()
        return

    # Workflow commands
    if not validate_workflow_arguments(cli_arguments):
        logger.info("Please check your command line arguments. Use 'spatialsnake --help' for more information")
        sys.exit(2)

    command = cli_arguments.get("<command>")
    if command in ['single_analysis', 'compare_analysis']:
        runner = WorkflowRunner(cli_arguments)
        runner.execute()
        return

    logger.error("No valid command was selected. Use 'spatialsnake --help' for more information")
    sys.exit(2)

if __name__ == '__main__':
    main()
