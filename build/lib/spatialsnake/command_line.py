#!/usr/bin/env python
'''
Created on 2025/10/5
spatialsnake main
@author: lzh,xulabgdpu,1714074171@qq.com
'''
import re
import warnings
warnings.filterwarnings("ignore")
from docopt import docopt
import os
import sys
import subprocess
import shutil
import datetime
import random
import timeit
import errno
import yaml
from yaml.loader import SafeLoader
from subprocess import call
import pathlib
from collections import defaultdict
import spatialsnake



spatialsnake_path=os.path.dirname(spatialsnake.__file__)
option = ["integrate","preprocess","clustering","annotion_help","annotion","compare_analyze","advance_analysis","splitting","merge","transform"]

__author__ = 'lzh'
__version__= '0.1.0'
__logo__="""
     
   _____
  /     \\    SpatialSnake
 |  ()  |   ------------
  \\  ^  /   Automated
   |||||     Spatial
  /|||||\\   Analysis
 |/|||||\\|  Pipeline
   ~~~~~    v0.1.0                                              
"""  

__licence__="""
MIT License
Copyright (c) 2025
...
"""

##### [options]  : 以--开头的所有参数通配符

__doc__=f"""Main spatialsnake executable, version: {__version__}
{__logo__} 

Usage:
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
    --runpipe <TEXT>          Run  which analysis analysis.[default: advance_analysis]
    --senic_input <DIR>       Input file for PySCENIC analysis.[default: sample.zarr]
    --motifs_input <FILE>     PySCENIC database directory.[default: motifs-v9-nr.hgnc-m0.001-o0.0.tbl]
    --feather_input <FILE>    path for necessary file of pyscenic.[default: hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather]
    --tfs_input <FILE>        path for necessary file of pyscenic.[default: hs_hgnc_tfs.txt]
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

"""

# Preprocessing Step Options (--option preprocess):
#     --integration-method TEXT   Integration method [default: harmony].
#     --annotion_list FILE    for the filter params in different sample
#     --min_cells INT         Minimum spots per gene [default: 3].
#     --min_genes INT         Minimum genes per spot [default: 200].
#     --variable BOOL         Filter the variable spot to analysis [default: False].
#     --harmony BOOL          harmony method [default: True].
#     --seg_filter BOOL       to seg filter the differnet sample dataset when command compare_anaysis.
#     --NEIGHBORS FLOAT       neighbors for pca umap.
# Clustering Step Options (--option clustering):
#     --resolution FLOAT   Cluster resolution [default: 0.5].
#     --cluster_algorithm TEXT Clustering algorithm [default: leiden].
#     --tsene BOOL        umap [default:False]
#     --MIN_DIST FLOAT    umap_key [default:0.3]
#     --SPREAD FLOAT      umap_key [default:1]
# 
# Annotation Help Step Options (--option annotion_help):
#     --image_slice BOOL        containing marker genes for cell types[default: False].
#     --markers_algorithm TEXT       Automatically detect marker genes [default: wilcoxon].
#     --shape_type TEXT         Automatically detect marker genes [default: cell_boundaries].
#     --image_type TEXT         Automatically detect marker genes [default: hires].
#     --spacies TEXT            Automatically detect marker genes [default: human].
#     --image_slice BOOL              params for the image slice to depandent size [default: False].
#     --x1 INT
#     --x2 INT
#     --y1 INT
#     --y2 INT
# Compare_analyze option Options (--option compare_analysis)    
#     --cell_focus TEXT         celltype you focus to compare in different sample.
#     --compare_algorithm TEXT  compare analysys [default: DEseq2].
# Annotation option Options (--option annotion):
#     --annotation-file FILE    Annotation file for cell typing (required for annotion step)
#     --anno_algorithm TEXT     Annotation method [default: mannul].
#     --shape_type TEXT         Automatically detect marker genes [default: cell_boundaries].
#     --image_type TEXT         Automatically detect marker genes [default: hires].
#     --slice BOOL              params for the image slice to depandent size [default: False].
#     --x1 INT
#     --x2 INT
#     --y1 INT
#     --y2 INT
#     --max_epochs_reference INT    params for cell2Location model train and test [default: 250].
#     --remove_mt BOOL              params for cell2Location model train and test [default: True].
#     --N_cells_per_location INT    params for cell2Location model train and test [default: 30].
#     --max_epochs_st INT           params for cell2Location model train and test [default: 30000].
#     --device TEXT                 cpu or GPU accelerate [default: cuda].
#     
# Advanced Analysis option Options (--option advance_analysis):
#     --advance_channel TEXT        Run  which analysis analysis.
#     --pyscenic-input FILE   Input file for PySCENIC analysis.
#     --pyscenic-db FILE      PySCENIC database directory.
#     --pyscenic-feature FILE path for necessary file of pyscenic.
#     --pyscenic-tfs FILE     path for necessary file of pyscenic.
#     --cell_attr TEXT        cell_id for pyscenic [default: cell_id]
#     --workers INT           workers for pyscenic [default: 8].
#     --count-data TEXT       gene type for cellPhoneDB [default: hgnc_symbol].
#     --threads INT           workers for cellphoneDB [default: 8].
#     --output_name           output name for cellPhoneDB [default: Normal].
#     
# Cell Segmentation Options (applicable to multiple option):
#     --zarr_file FILE        seg with the sample name or region

def check_command_line_arguments(arguments):
    if arguments["useful_tool"]:
        return True
    if not os.path.exists(arguments["<INPUT>"]):
        print("sample list file not found ",arguments["<INPUT>"])
        return False
    if arguments["--option"] in ["integrate","clustering","annotion_help","compare_analyze","advance_analysis"] and arguments["<INPUT>"]!="sample.txt":
        print('.please confirm the file name are sample_list.txt or annotion_list or filter_list')
        return False
    if arguments["<TYPE>"] not in ['visium','visium_segment','visium_HD','xenium','Merfish','slide_seq']:
        print("please select the correct spatialdata type like:'visium','visium_segment','visium_HD','xenium','Merfish','slide_seq'")
        return False
    if arguments["--option"] not in option:
        print("your option are not correct please select the correct step to analysis or not select the option to run the minimize step of analysis")
        print("correct option include:"+option)
        return False
    if "--configfile" in arguments and arguments["--configfile"]!='config.yaml':
        if not  os.path.isfile(arguments["--configfile"]):
          print("please select a .yaml file for --configfile")
          return False
    return True



class CommandLine:
    def __init__(self):
        self.snakemake="snakemake --rerun-incomplete -k "
        self.runid="".join(random.choices("abcdefghisz",k=3) + random.choices("123456789",k=5))
        self.config=[]
        self.configfile_loaded=False
        self.is_useful_tool_sample=False
        self.is_this_an_useful_tool_run=False
        self.parameters=dict()
        self.log=True
        
    def __str__(self):
        return self.snakemake
    def __repr__(self):
        return self.snakemake
    

    def add_config_argument(self):
        self.snakemake = self.snakemake + " --config " + " ".join(self.config)


    def load_configfile_if_available(self,arguments):
        step = arguments["--option"] if not self.is_this_an_useful_tool_run else arguments["<command>"]
        if self.configfile_loaded is False:
            if "--configfile" in arguments and os.path.isfile(arguments["--configfile"]):
                self.snakemake = self.snakemake + " --configfile={}".format(arguments["--configfile"])
                configfile=arguments["--configfile"]
                self.configfile_loaded=True
            else:
               self.snakemake = self.snakemake + " --configfile={}".format(spatialsnake_path + f"/workflow/envs/{step}.yaml")
               configfile=spatialsnake_path + f"/workflow/envs/{step}.yaml"
               arguments["--configfile"]=spatialsnake_path + f"/workflow/envs/{step}.yaml"
            with open(configfile) as f:
               self.parameters=yaml.load(f,Loader=SafeLoader)

    def prepare_arguments(self,arguments):
        jobs = arguments.get('--jobs', 4)
        self.snakemake = self.snakemake + " -j {} ".format(jobs)
        self.snakemake = self.snakemake +  " -s {} ".format(f"{spatialsnake_path}/workflow/Snakefile")
        self.load_configfile_if_available(arguments)
        if arguments['--option'] in ["integrate","preprocess","clustering","annotion_help","compare_analyze","all"]:
          self.config.append("sample_list={}".format(arguments['<INPUT>']))
        self.config.append(f"spatialsnake_path={spatialsnake_path}/")
        for i,b in arguments.items():
            if i not in ["--jobs","--configfile","--option","--unlock","--remove","--dry","--help","--version","<INPUT>","<command>","--install-packages","<TYPE>","useful_tool","<INTEGRATED_FILE>"]:
                k=i.lstrip("--")
                if k in ["min_cells", "min_genes", "x1", "x2", "y1", "y2", "workers", "threads"]:
                  try:
                    int(b)
                  except ValueError:
                    print(f"Error: {k} must be an integer")
                    sys.exit(1)
                if self.parameters.get(k)==None:
                  continue
                if self.configfile_loaded is False: 
                    self.config.append(k + "=" + str(b))
                    self.parameters[k]=str(b)
                else:
                    if self.parameters.get(k) and i not in sys.argv:
                        self.config.append(k + "=" + str(self.parameters.get(k)))
                    else:
                        self.config.append(k + "=" + str(b))
                        self.parameters[k]=str(b)
        
        self.config.append("runid={}".format(self.runid))
        if arguments["--option"]:
            self.config.append("option={}".format(arguments["--option"]))
        if self.is_this_an_useful_tool_run is False:
            self.config.append("channel={}".format(arguments['<command>']))
            self.config.append("run_type={}".format(arguments["<TYPE>"]))
        if "--dry" in arguments and arguments["--dry"]:
            self.snakemake = self.snakemake + " -n "
            self.log=False
        if "--unlock" in arguments and arguments["--unlock"]:
            self.snakemake = self.snakemake + " --unlock "
            self.log=False
        if "--remove" in arguments and arguments["--remove"]:
            self.snakemake = self.snakemake + " --delete-all-output "
            self.log=False
        self.add_config_argument()
        
    
    def write_to_log(self,start,arguments):
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        logname = f"spatialsnake_{self.runid}_{timestamp}_runlog.log"
        stop = timeit.default_timer()
        if self.log:
            with open(logname,"w") as f:
                f.write(__logo__ + "\n")
                f.write("Run ID : " + self.runid + "\n")
                f.write("spatialsnake version : " + __version__ + "\n")
                f.write("spatialsnake arguments : " + " ".join(sys.argv) + "\n")
                f.write("the running step : "+self.parameters["option"]+"\n")
                f.write("------------------------------" + "\n")
                f.write("Snakemake arguments : " + str(self.snakemake) + "\n\n")
                f.write("------------------------------" + "\n")
                f.write("Run parameters in this option:\n")
                for key, value in sorted(self.parameters.items()):
                    if value is not None and value != "":
                        f.write(f"  {key.ljust(25)} {value}\n")
                f.write("\n")
                f.write("Total run time: {t:.2f} mins \n".format(t=(stop-start)/60))
                f.write("Useful Information:\n")
                f.write("-" * 20 + "\n")
                f.write("  ⚠  For help: [spatialsnake --help ]\n")
                f.write("  ⚠  For setting more params please run: [spatialsnake produce-file --option=step]\n")
                f.write("  ⚠  Output files are stored in the 'results' directory\n")
                f.write("  ⚠  Use 【spatialsnake [command] [..] --config-file config.yaml】 to run spatialsnake if you want to customize the parameter settings in config.yaml file. \n")
                f.write("=" * 60 + "\n")
                
class CommandLine_useful_tools:
    def __init__(self):
        self.snakemake="python "
        self.runid="".join(random.choices("abcdefghisz",k=3) + random.choices("123456789",k=5))
        self.config=[]
        self.configfile_loaded=False
        self.is_this_an_useful_tool_run=True
        self.parameters=dict()
        self.log=True
        
    def __str__(self):
        return self.snakemake
    def __repr__(self):
        return self.snakemake
      
    def add_config_argument(self):
        self.snakemake = self.snakemake + " ".join(self.config)

    def load_configfile_if_available(self,arguments):
        if arguments["--option"]:
            tool = arguments["--option"]
            # self.config.append("option={}".format(arguments["--option"]))
            self.config.append(spatialsnake_path + f"/workflow/function/{tool}.py")
        print(spatialsnake_path + f"/workflow/envs/{tool}.yaml")
        if self.configfile_loaded is False:
            print(spatialsnake_path)
            if "--configfile" in arguments and os.path.isfile(arguments["--configfile"]):
                configfile=arguments["--configfile"]
                self.configfile_loaded=True
                print(configfile)
            else:
               configfile=spatialsnake_path + f"/workflow/envs/{tool}.yaml"
               arguments["--configfile"]=spatialsnake_path + f"/workflow/envs/{tool}.yaml"
               print(configfile)
            with open(configfile) as f:
               self.parameters=yaml.load(f,Loader=SafeLoader)
    def prepare_arguments(self,arguments):
        # jobs = arguments.get('--jobs', 4)
        self.load_configfile_if_available(arguments)
        self.config.append("--INPUT {}".format(arguments['<INPUT>']))
        for i,b in arguments.items():
            if i not in ["--jobs","--configfile","--option","useful_tool","<INPUT>"]:
                k=i.lstrip("--")
                if self.parameters.get(k)==None:
                  continue
                if self.configfile_loaded is False: 
                    self.config.append(i + " " + str(b))
                    self.parameters[k]=str(b)
                else:
                    if self.parameters.get(k) and i not in sys.argv:
                        self.config.append(i +" "+ str(self.parameters.get(k)))
                    else:
                        self.config.append(i + " "+str(b))
                        self.parameters[k]=str(b)
        self.add_config_argument()
        
    
    def write_to_log(self,start,arguments):
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        logname = f"spatialsnake_{self.runid}_{timestamp}_runlog.log"
        stop = timeit.default_timer()
        if self.log:
            with open(logname,"w") as f:
                f.write(__logo__ + "\n")
                f.write("Run ID : " + self.runid + "\n")
                f.write("spatialsnake version : " + __version__ + "\n")
                f.write("spatialsnake arguments : " + " ".join(sys.argv) + "\n")
                f.write("the running step : "+arguments["--option"]+"\n")
                f.write("------------------------------" + "\n")
                f.write("Snakemake arguments : " + str(self.snakemake) + "\n\n")
                f.write("------------------------------" + "\n")
                f.write("Run parameters in this option:\n")
                for key, value in sorted(self.parameters.items()):
                    if value is not None and value != "":
                        f.write(f"  {key.ljust(25)} {value}\n")
                f.write("\n")
                f.write("Total run time: {t:.2f} mins \n".format(t=(stop-start)/60))
                f.write("Useful Information:\n")
                f.write("-" * 20 + "\n")
                f.write("  ⚠  For help: [spatialsnake --help ]\n")
                f.write("  ⚠  For setting more params please run: [spatialsnake produce-file --option=step]\n")
                f.write("  ⚠  Output files are stored in the 'results' directory\n")
                f.write("  ⚠  Use 【spatialsnake [command] [..] --config-file config.yaml】 to run spatialsnake if you want to customize the parameter settings in config.yaml file. \n")
                f.write("=" * 60 + "\n")



 
def run_useful_tool(arguments):
    start = timeit.default_timer()
    snakemake_argument=CommandLine_useful_tools()
    snakemake_argument.is_this_an_useful_tool_run = True
    snakemake_argument.prepare_arguments(arguments)
    print(snakemake_argument)
    subprocess.check_call(str(snakemake_argument),shell=True)
    snakemake_argument.write_to_log(start,arguments)

def run_workflow(arguments):
    start = timeit.default_timer()
    snakemake_argument=CommandLine()
    snakemake_argument.prepare_arguments(arguments)
    subprocess.check_call(str(snakemake_argument),shell=True)
    snakemake_argument.write_to_log(start,arguments)


def main():
        cli_arguments = docopt(__doc__, version=__version__)
        if cli_arguments["produce-file"]:
            step = cli_arguments["--option"] if not cli_arguments['<command>'] == 'useful_tool' or not cli_arguments['<command>'] == 'transform' else cli_arguments["<command>"]
            if step not in ["integrate","preprocess","clustering","annotion_help","annotion","compare_analyze","advance_analysis","all"]:
              print("please setting correct params : --option=<step_name>   or  --option=all to get all step params")
              return
            print(f"Generating config.yaml file: {step}.yaml..........")
            print("You can use this as a config-file for a spatialsnake run. You may change the settings in it.")
            if step=="all":
              shutil.copyfile(spatialsnake_path + "/config.yaml", 'config.yaml')
            else:
              shutil.copyfile(spatialsnake_path + f"/workflow/envs/{step}.yaml", f'{step}.yaml')
            print("⚠     How to setting your own params:")
            print("Add params: --config-file <file-path> in the command line when you run the pipeline")
            return
        elif cli_arguments["install-packages"]:
            r_script_path = spatialsnake_path + "/workflow/scripts/install_packages.R"
            subprocess.check_call(["Rscript", r_script_path])
            return
        if not check_command_line_arguments(cli_arguments):
            print("""Please check your command line arguments. Use "spatialsnake --help" for more information""")
            return
        if cli_arguments['<command>'] == 'single_analysis':
            run_workflow(cli_arguments)
        elif cli_arguments['<command>'] == 'compare_analysis':
            run_workflow(cli_arguments)
        elif cli_arguments['useful_tool']:
            run_useful_tool(cli_arguments)
        elif cli_arguments['transform']:
           run_transform(cli_arguments)
