#!/usr/bin/env python
'''
spatialsnake main
@author: Zhenghao Lin,2400507123@gdpu.edn.stu
'''
import datetime
import importlib.util
import logging
import os
import random
import shutil
import subprocess
import sys
import timeit
import warnings
from typing import Any, Dict, List

import spatialsnake
import yaml
from docopt import docopt
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

# Constants
SPATIALSNAKE_PATH = os.path.dirname(spatialsnake.__file__)
VALID_OPTIONS = [
    "integrate", "preprocess", "clustering", "annotation_help", "annotation",
    "compare_stage", "advance_analysis", "reclustering", "splitting", "merge", "transform"
]

__author__ = 'lzh'
__version__ = '0.0.1'
__logo__ = """

  ╭─── SpatialSnake · v0.2.4 ───╮
  │                              │
  │    ●───●───●───●───●───●    │
  │    │ ╲ │ ╱ │ ╲ │ ╱ │ ╲ │    │
  │    ●───●───●───●───●───●    │
  │                              │
  │   Spatial Transcriptomics    │
  │      Analysis Pipeline       │
  ╰──────────────────────────────╯
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
    reclustering
    annotation_help
    annotation
    compare_stage
    advance_analysis

Type Arguments:
    visium
    visium_segment
    visium_HD
    xenium
    Merfish
    slide_seq
    stereoseq

INPUT Arguments:
    sample.txt
    annotation.txt
    filter_list
    
Basic Configuration:
    --configfile <FILE>    Configuration file in YAML format [default: config.yaml].

Integration Step Options (--option integrate):
    --cells_boundaries <BOOL>    xenium key in load in data [default: False].
    --nucleus_boundaries <BOOL>  xenium key in load in data [default: False].
    --nucleus_labels <BOOL>      xenium key in load in data [default: False].
    --morphology_mip <BOOL>      xenium key in load in data [default: False].
    --bin_size <INT>             legacy fallback for Stereo-seq; prefer `sample.txt` column 3 `input_spec` (`cellbin`, `adjusted_cellbin`, or `50,150`) [default: 50].
    --merscope_z_layers <TEXT>   optional z layers for `spatialdata_io.merscope`, e.g. `0` or `0,1,2`.
    --merscope_region_name <TEXT> optional region name for `spatialdata_io.merscope`.
    --merscope_transcripts <BOOL> load transcripts in `spatialdata_io.merscope` [default: True].
    --merscope_cells_boundaries <BOOL> load cell boundaries in `spatialdata_io.merscope` [default: True].
    --merscope_cells_table <BOOL> load cells table in `spatialdata_io.merscope` [default: True].
    --merscope_mosaic_images <BOOL> load mosaic images in `spatialdata_io.merscope` [default: True].

Preprocessing Step Options (--option preprocess):
    --min_cells <INT>         Minimum spots per gene [default: 3].
    --min_genes <INT>         Minimum genes per spot [default: 200].
    --seg_filter <BOOL>       to seg filter the differnet sample dataset when command compare_anaysis [default: False].
    --filter_list <FILE>      filename of filter [default: False]
    --batch_method <TEXT>     batch method for multiple sample analysis [default: harmony]
    --sketch <BOOL>           whether use sketch method to analysis [default: False]
    --mt_threshold <FLOAT>    the mt params percent to filter the cell [default: 50.0].
    
Clustering Step Options (--option clustering):
    --resolution <FLOAT>        Cluster resolution [default: 0.5].
    --cluster_algorithm <TEXT>  Clustering algorithm [default: leiden].
    --tsene <BOOL>              umap [default:False].
    --n_clusters <INT>          kmeans params of cluster [default: 15].
    --pcs <INT>                 dimension pca select [default: 25].

Reclustering Step Options (--option reclustering):
    --recluster_resolution <FLOAT>      Leiden resolution [default: 0.8].
    --recluster_n_top_genes <INT>       Highly variable genes count [default: 2000].
    --recluster_neighbors <INT>         Neighbors for graph construction [default: 15].
    --recluster_n_pcs <INT>             Number of PCs for neighbors [default: 30].
    --recluster_marker_method <TEXT>    Marker test method [default: wilcoxon].
    --recluster_min_pct <FLOAT>         Min fraction for marker filtering [default: 0.1].
    --recluster_logfc_threshold <FLOAT> Min log2FC for marker filtering [default: 0.25].
    
Annotation Help Step Options (--option annotation_help):
    --markers_algorithm <TEXT>       Automatically detect marker genes [default: wilcoxon].
    --spacies <TEXT>            Automatically detect marker genes [default: human].
Compare_stage option Options (--option compare_stage)
    --cell_focus <TEXT>         celltype you focus to compare in different sample[default: None].
    --compare_algorithm <TEXT>  compare analysys [default: DEseq2].
Annotation option Options (--option annotation):
    --annotation-file <FILE>    Annotation file for cell typing (required for annotation step)
    --anno_algorithm <TEXT>     Annotation method (mannul/reannotation/cell2Location/RCTD) [default: mannul].
    --shape_type <TEXT>         Automatically detect marker genes [default: cell_boundaries].
    --image_type <TEXT>         Automatically detect marker genes [default: hires].
    --vis_mode <TEXT>           Spatial visualization mode [auto/point/shape, default: auto].
    --point_size <FLOAT>        Point render size for point-based visualization [default: 2.5].
    --device <TEXT>                 cpu or GPU accelerate [default: cuda].
    --max_cores <INT>               max cores for parallel [default: 16].
    --zarr_input <DIR>       Input file for PySCENIC analysis.

Advanced Analysis option Options (--option advance_analysis):
    --runpipe <TEXT>          Run  which analysis analysis.[default: advance_analysis]
    --senic_input <DIR>       Input file for PySCENIC analysis.[default: sample.zarr]
    --motifs_input <FILE>     PySCENIC database directory.
    --feather_input <FILE>    path for necessary file of pyscenic.
    --tfs_input <FILE>        path for necessary file of pyscenic.
    --count-data <TEXT>       gene type for cellPhoneDB [default: hgnc_symbol].
    --threads <INT>           workers for cellphoneDB [default: 16].
    --output_name <TEXT>      output name for cellPhoneDB [default: Normal].
    --workers <INT>           workers for pysenic [default: 32].
    --niche_col <TEXT>          niche column in the zarr/table/obs

useful_tool splitting Option:
    --output_dir=<TEXT>       output dir for splitted file [default: results/useful_results]
    --split_by=<TEXT>         slice out with the barcode in table[anndata] .obs [default: clusters]
    --barcodes=<TEXT>         comma-separated values for split_by filter [default: ""]
    --roi_csv=<TEXT>          csv file or directory for ROI splitting [default: ""]
    --shape_elements=<TEXT>   slice out with the shape [default: None]
    --max_x=<FLOAT>         coordinate of image boundaries [default: 0]
    --min_x=<FLOAT>         coordinate of image boundaries [default: 2000]
    --max_y=<FLOAT>         coordinate of image boundaries [default: 2000]
    --min_y=<FLOAT>         coordinate of image boundaries [default: 0]

useful_tool merge Option:
    --merge_by=<TEXT>               merge by cluster celltype sample or reannotation [default: sample]
    --reordering=<BOOL>             whether reordering the cluster when concat the [default: False]
    --re_sample=<BOOL>              whether add the sample lable[default: False]
    --cluster_key=<TEXT>            the concat lable in the zarr/table/obs [default: clusters]
    --annotation_csv=<TEXT>         csv path, csv directory, or comma-separated csv paths [default: ""]
    --csv_cell_col=<TEXT>           cell id column in annotation csv [default: Barcode]
    --csv_label_col=<TEXT>          label column in annotation csv [default: Grouped_Annotation]
    --input_cell_col=<TEXT>         cell id column in base zarr table obs [default: cell_id]
    --target_col=<TEXT>             output column name written to base zarr [default: sub_celltype]
    --original_celltype_col=<TEXT>  original celltype column used when first creating target_col [default: celltype]

useful_tool transform Option:
    --save_image=<BOOL>            save images in h5ad [default: True]
    --transform_from=<TEXT>        transform from [default: zarr]
    --transform_to=<TEXT>          transform to [default: h5ad]

General Options:
    -j <INT>, --jobs <INT>   Number of CPU cores [default: 16].
    --results_folder <DIR>     Output directory [default: results].

Utility Options:
    --install-packages   Install required packages.
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

    def write_to_log(self, start_time: float):
        """Write execution details to log file."""
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_dir = "log"
        logname = os.path.join(log_dir, f"spatialsnake_{self.runid}_{timestamp}_runlog.log")
        stop_time = timeit.default_timer()

        if self.log_enabled:
            try:
                os.makedirs(log_dir, exist_ok=True)
                with open(logname, "w") as f:
                    f.write(__logo__ + "\n")
                    f.write(f"Run ID : {self.runid}\n")
                    f.write(f"spatialsnake version : {__version__}\n")
                    f.write(f"spatialsnake arguments : {' '.join(sys.argv)}\n")
                    option_val = self.arguments.get("--option", "unknown")
                    f.write(f"the running step : {option_val}\n")
                    f.write("-" * 30 + "\n")
                    f.write(f"Command arguments : {self.cmd_str}\n\n")
                    f.write("-" * 30 + "\n")
                    f.write("Run parameters in this option:\n")
                    for key, value in sorted(self.parameters.items()):
                        if value is not None and value != "":
                            f.write(f"  {key.ljust(25)} {value}\n")
                    f.write("\n")
                    f.write(f"Total run time: {(stop_time - start_time) / 60:.2f} mins \n")
                    f.write("Useful Information:\n")
                    f.write("-" * 20 + "\n")
                    f.write("  ⚠  For help: [spatialsnake --help ]\n")
                    f.write("  ⚠  For setting more params please run: [spatialsnake produce-file --option=step]\n")
                    f.write("  ⚠  Output files are stored in the 'results' directory\n")
                    f.write("  ⚠  Use 【spatialsnake [command] [..] --config-file config.yaml】 to run spatialsnake if you want to customize the parameter settings in config.yaml file. \n")
                    f.write("=" * 60 + "\n")
            except Exception as e:
                logger.error(f"Failed to write log file: {e}")

    def execute(self):
        """Execute the prepared command."""
        start_time = timeit.default_timer()
        try:
            self.prepare_arguments()
            logger.info(f"Executing command: {self.cmd_str}")
            subprocess.check_call(str(self.cmd_str), shell=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Command execution failed with return code {e.returncode}")
            sys.exit(e.returncode)
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            sys.exit(1)
        finally:
            self.write_to_log(start_time)


class WorkflowRunner(BaseCommandLine):
    """Handles Snakemake workflow execution."""

    def __init__(self, arguments: Dict[str, Any]):
        super().__init__(arguments)
        self.snakemake_cmd = "snakemake --rerun-incomplete -k "
        self.cmd_str = self.snakemake_cmd

    def add_config_argument(self):
        self.cmd_str += " --config " + " ".join(self.config)

    def load_configfile(self):
        step = self.arguments.get("--option")
        configfile = None
        if self.arguments.get("--configfile") and os.path.isfile(self.arguments["--configfile"]):
            configfile = self.arguments["--configfile"]
            self.cmd_str += f" --configfile={configfile}"
            self.configfile_loaded = True
        else:
            # Fallback to default env config
            default_config = os.path.join(SPATIALSNAKE_PATH, f"workflow/envs/{step}.yaml")
            self.cmd_str += f" --configfile={default_config}"
            configfile = default_config
            self.arguments["--configfile"] = default_config
            
        if configfile and os.path.exists(configfile):
            try:
                with open(configfile) as f:
                    self.parameters = yaml.load(f, Loader=SafeLoader) or {}
            except Exception as e:
                logger.warning(f"Failed to load config file {configfile}: {e}")

    def prepare_arguments(self):
        jobs = self.arguments.get('--jobs', 4)
        self.cmd_str += f" -j {jobs} "
        self.cmd_str += f" -s {os.path.join(SPATIALSNAKE_PATH, 'workflow/Snakefile')} "
        
        self.load_configfile()
        
        if self.arguments.get('--option') in ["integrate", "preprocess", "clustering", "reclustering", "annotation_help", "annotation", "compare_stage", "all", "advance_analysis"]:
            self.config.append(f"sample_list={self.arguments.get('<INPUT_FILE>')}")
            
        self.config.append(f"spatialsnake_path={SPATIALSNAKE_PATH}/")
        
        # Parse dynamic arguments
        exclude_keys = {
            "--jobs", "--configfile", "--option", "--unlock", "--remove", 
            "--dry", "--help", "--version", "<INPUT_FILE>", "<command>", 
            "--install-packages", "<TYPE>", "useful_tool", "<INTEGRATED_FILE>"
        }
        for key, value in self.arguments.items():
            if key in exclude_keys:
                continue
                
            clean_key = key.lstrip("--")
            if clean_key == "annotation-file":
                clean_key = "annotation_list"
            
            # Integer validation
            if clean_key in ["min_cells", "min_genes", "x1", "x2", "y1", "y2", "workers", "threads", "bin_size"]:
                try:
                    int(value)
                except (ValueError, TypeError):
                    logger.error(f"Error: {clean_key} must be an integer, got '{value}'")
                    sys.exit(1)
            
            # Check if parameter is relevant for loaded config
            if self.parameters.get(clean_key) is None:
                continue

            if not self.configfile_loaded:
                if value is None:  ### 测试专用
                    continue
                self.config.append(f"{clean_key}={value}")
                self.parameters[clean_key] = str(value)
            else:
                # Prefer value from config file unless overridden in sys.argv?
                # The original logic was a bit weird: 
                # if self.parameters.get(k) and i not in sys.argv: use config value
                # else: use argument value
                # But 'i' is the key from arguments dict (e.g. "--min_cells").
                
                if self.parameters.get(clean_key) and key not in sys.argv:
                     self.config.append(f"{clean_key}={self.parameters[clean_key]}")
                else:
                     self.config.append(f"{clean_key}={value}")
                     self.parameters[clean_key] = str(value)
        self.config.append(f"runid={self.runid}")
        
        if self.arguments.get("--option"):
            self.config.append(f"option={self.arguments['--option']}")
            
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
        if self.arguments.get("--configfile") and os.path.isfile(self.arguments["--configfile"]):
            configfile = self.arguments["--configfile"]
            self.configfile_loaded = True
            logger.info(f"Using config file: {configfile}")
        else:
            default_config = os.path.join(SPATIALSNAKE_PATH, f"workflow/envs/{tool}.yaml")
            configfile = default_config
            self.arguments["--configfile"] = default_config
            logger.info(f"Using default config file: {configfile}")
            
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
                self.config.extend(inputs)
            else:
                 self.config.append(str(inputs))
        elif self.arguments.get("--option") == "splitting":
            inputs = self.arguments.get('<INPUT>')
            input_value = inputs[0] if isinstance(inputs, list) and len(inputs) > 0 else inputs
            self.config.append(f"--INPUT_FIlE {input_value}")
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
        
        exclude_keys = {"--jobs", "--configfile", "--option", "useful_tool", "<INPUT>"}
        
        for key, value in self.arguments.items():
            if key in exclude_keys:
                continue
                
            clean_key = key.lstrip("--")
            if self.parameters.get(clean_key) is None:
                continue
                
            if not self.configfile_loaded:
                self.config.append(f"{key} {value}")
                self.parameters[clean_key] = str(value)
            else:
                if self.parameters.get(clean_key) and key not in sys.argv:
                    self.config.append(f"{key} {self.parameters[clean_key]}")
                else:
                    self.config.append(f"{key} {value}")
                    self.parameters[clean_key] = str(value)
                    
        self.add_config_argument()


def validate_workflow_arguments(arguments: Dict[str, Any]) -> bool:
    """Validate arguments for workflow commands."""
    input_file = arguments.get("<INPUT_FILE>")
    if input_file and not os.path.exists(input_file):
        logger.error(f"Sample list file not found: {input_file}")
        return False
        
    option = arguments.get("--option")
    if option in ["integrate", "clustering", "reclustering", "annotation_help", "compare_stage", "advance_analysis"]:
        if input_file != "sample.txt":
             # This seems like a strict requirement in the original code
             logger.warning("Please confirm the file name is sample.txt (or appropriate list file)")
             # The original code returned False here if not "sample.txt"
             # "if ... and arguments["<INPUT_FILE>"]!="sample.txt": return False"
             # I will keep it strict as per original.
             logger.error("For integrate/clustering/etc., input file must be 'sample.txt'")
             return False

    type_arg = arguments.get("<TYPE>")
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
    valid_types = ['visium', 'visium_segment', 'visium_HD', 'xenium', 'Merfish', 'slide_seq', 'stereoseq']
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

    if option == "splitting":
        # Original code used arguments["split_by"] but docopt usually keeps --. 
        # Wait, the usage says `[--split_by=<TEXT>]`. Docopt key would be `--split_by`.
        # But original code accessed `arguments["split_by"]`. This implies some preprocessing or docopt implementation detail?
        # Standard docopt returns keys as specified in Usage or Options.
        # In Usage: `spatialsnake useful_tool ...`
        # In Options: `--split_by=<TEXT>`
        # So key is `--split_by`.
        # However, `CommandLine.prepare_arguments` iterates `arguments.items()` and strips `--`.
        # `check_arguments_inputfile` accessed `arguments["split_by"]`.
        # If the key is `--split_by`, `arguments["split_by"]` would raise KeyError.
        # Unless the user passed `split_by` as a command/argument? No.
        # I suspect the original code might have had issues or I am missing something about docopt.
        # I will assume the key is `--split_by`.
        
        val = arguments.get("--split_by")
        # Also need to check if user passed it?
        # If user didn't pass it, docopt uses default "clusters".
        if val and val not in ["sample", "image", "cluster", "clusters","group","region","celltype","ROI","ROIs"]:
             # "cluster" was in original check list: ["sample",'image',"cluster"]
             # But default is "clusters".
             logger.error("split_by must be one of: sample, image, cluster, ROI")
             return False

    return True


def install_packages():
    _install_pybanksy_if_needed()
    """Install R and optional Python packages."""
    r_script_path = os.path.join(SPATIALSNAKE_PATH, "workflow/scripts/install_packages.R")
    try:
        logger.info(f"Running R script: {r_script_path}")
        subprocess.check_call(["Rscript", r_script_path])
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install R packages: {e}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Rscript not found or script missing.")
        sys.exit(1)


def _install_pybanksy_if_needed():
    """pybanksy metadata declares numpy<2.0 which conflicts with our numpy>=2.
    Install it with --no-deps since all its true dependencies are already
    satisfied by spatialsnake core requirements."""
    try:
        if importlib.util.find_spec("pybanksy") is not None:
            return
    except Exception:
        return

    logger.info("Installing pybanksy (BANKSY spatial clustering) ...")
    try:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "pybanksy==1.3.4", "--no-deps"],
        )
        logger.info("pybanksy installed successfully.")
    except subprocess.CalledProcessError:
        logger.warning(
            "Failed to install pybanksy automatically. "
            "You can install it manually later: pip install pybanksy==1.3.4 --no-deps"
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
        return

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
        logger.info("⚠  How to set your own params:")
        logger.info("   Add params: --config-file <file-path> in the command line")
        
    except FileNotFoundError:
        logger.error(f"Source config file not found: {src}")
    except Exception as e:
        logger.error(f"Error generating config file: {e}")


def main():
    """Main entry point."""
    try:
        cli_arguments = docopt(__doc__, version=__version__)
    except Exception as e:
        logger.error(f"Error parsing arguments: {e}")
        return

    # Debug print as in original? Maybe remove for production, but user asked to keep core functions.
    # Original code had `print(cli_arguments)` in some branches.
    # I'll rely on logging.

    if cli_arguments.get("produce-file"):
        generate_config_file(cli_arguments)
        return

    if cli_arguments.get("install-packages"):
        install_packages()
        return

    if cli_arguments.get("useful_tool"):
        # print(cli_arguments) # Original had this
        if not validate_tool_arguments(cli_arguments):
             return
        runner = ToolRunner(cli_arguments)
        runner.execute()
        return

    # Workflow commands
    # print(cli_arguments) # Original had this
    
    if not validate_workflow_arguments(cli_arguments):
        logger.info("Please check your command line arguments. Use 'spatialsnake --help' for more information")
        return

    command = cli_arguments.get("<command>")
    if command in ['single_analysis', 'compare_analysis']:
        runner = WorkflowRunner(cli_arguments)
        runner.execute()

if __name__ == '__main__':
    main()
