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
option = ["integrate","preprocess","clustering","annotion_help","annotion","compare_analyze","advance_analysis"]

__author__ = 'lzh'
__version__= '0.1.0'
__logo__="""
     
   _____
  /     \\    SpatialSnake
 |  ()  |   ------------
  \\  ^  /   Automated
   |||||    Spatial
  /|||||\\   Analysis
 |/|||||\\|  Pipeline
   ~~~~~    v0.1.0                                              
"""  

__licence__="""
MIT License
Copyright (c) 2025 lzh  1714074171@qq.com
...
"""



__doc__=f"""Main spatialsnake executable, version: {__version__}
{__logo__} 

Usage:
    spatialsnake <command> <INPUT> <TYPE> [--option=options] [params]
    spatialsnake segout <INTEGRATED_ZARR> [--option=options] [params]
    spatialsnake produce-file
    spatialsnake install-packages
    spatialsnake (-h | --help)
    spatialsnake --version

Main Commands:
    single_analysis      Process single spatial transcriptomics dataset (runs all basic steps except advance_analysis by default)
    compare_analysis     Compare multiple spatial transcriptomics datasets
    segout               Split integrated ZARR data back into individual samples

option Control:
    --option OPTION          Run specific analysis step (for single_analysis and compare_analysis)
                         Options: integrate, preprocess, clustering, annotion_help, annotion, advance_analysis
                         Default: run all basic steps (integrate→preprocess→clustering→annotion_help→annotion)

Input Arguments:
    sample_list          Text file containing sample paths (one per line)
    segout_zarr      Integrated ZARR file to split (for segout command)
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

Integration Step Options (--option integrate):
    --SAMPLE_LIST FILE        sample list
    --integration-method TEXT   Integration method [default: harmony].
    --cells_boundaries BOOL    xenium key in load in data [default: False].
    --nucleus_boundaries BOOL  xenium key in load in data [default: False].
    --nucleus_labels BOOL      xenium key in load in data [default: False].
    --morphology_mip BOOL      xenium key in load in data [default: False].

Preprocessing Step Options (--option preprocess):
    --annotion_list FILE    for the filter params in different sample
    --min_cells INT         Minimum spots per gene [default: 3].
    --min_genes INT         Minimum genes per spot [default: 200].
    --variable BOOL         Filter the variable spot to analysis [default: False].
    --harmony BOOL          harmony method [default: True].
    --seg_filter BOOL       to seg filter the differnet sample dataset when command compare_anaysis.
    --NEIGHBORS FLOAT       neighbors for pca umap.
Clustering Step Options (--option clustering):
    --resolution FLOAT   Cluster resolution [default: 0.5].
    --cluster_algorithm TEXT Clustering algorithm [default: leiden].
    --tsene BOOL        umap [default:False]
    --MIN_DIST FLOAT    umap_key [default:0.3]
    --SPREAD FLOAT      umap_key [default:1]

Annotation Help Step Options (--option annotion_help):
    --image_slice BOOL        containing marker genes for cell types[default: False].
    --markers_algorithm TEXT       Automatically detect marker genes [default: wilcoxon].
    --shape_type TEXT         Automatically detect marker genes [default: cell_boundaries].
    --image_type TEXT         Automatically detect marker genes [default: hires].
    --spacies TEXT            Automatically detect marker genes [default: human].
    --image_slice BOOL              params for the image slice to depandent size [default: False].
    --x1 INT
    --x2 INT
    --y1 INT
    --y2 INT
Compare_analyze option Options (--option compare_analysis)    
    --cell_focus TEXT         celltype you focus to compare in different sample.
    --compare_algorithm TEXT  compare analysys [default: DEseq2].
Annotation option Options (--option annotion):
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
    
Advanced Analysis option Options (--option advance_analysis):
    --advance_channel TEXT        Run  which analysis analysis.
    --pyscenic-input FILE   Input file for PySCENIC analysis.
    --pyscenic-db FILE      PySCENIC database directory.
    --pyscenic-feature FILE path for necessary file of pyscenic.
    --pyscenic-tfs FILE     path for necessary file of pyscenic.
    --cell_attr TEXT        cell_id for pyscenic [default: cell_id]
    --workers INT           workers for pyscenic [default: 8].
    --count-data TEXT       gene type for cellPhoneDB [default: hgnc_symbol].
    --threads INT           workers for cellphoneDB [default: 8].
    --output_name           output name for cellPhoneDB [default: Normal].
    
Cell Segmentation Options (applicable to multiple option):
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

"""


def check_command_line_arguments(arguments):
    if not os.path.exists(arguments["<INPUT>"]):
        print("sample list file not found ",arguments["<INPUT>"])
        return False
    if arguments["--option"] in ["integrate","clustering","annotion_help","compare_analyze","advance_analysis"] and arguments["<INPUT>"]!="sample.txt":
        print('.please confirm the file name are sample_list.txt or annotion_list or filter_list')
        return False
    if arguments["segout"] and not os.path.isdir(arguments["segout_zarr"]):
        print("please select a zarr file to seg out")
        return False
    if arguments["<TYPE>"] not in ['visium','visium_segment','visium_HD','xenium','Merfish','slide_seq','Spatial_ATAC']:
        print("please select the correct spatialdata type like:'visium','visium_segment','visium_HD','xenium','Merfish','slide_seq','Spatial_ATAC'")
        return False
    if arguments["--option"] not in option:
        print("your option are not correct please select the correct step to analysis or not select the option to run the minimize step of analysis")
        print("correct option include:"+option)
        return False
    if arguments["--option"] and os.path.isdir(arguments["<INPUT>"]):
        print("You are running integrated option but you provided a directory, not a Seurat object file !")
        print("The default Seurat object is usually here, analyses_integrated/seurat/integrated.rds")
        return False
    if "--configfile" in arguments and arguments["--configfile"]:
         if not os.path.isfile(arguments["--configfile"]):
            print("Config file given not found : ",arguments["--configfile"])
            return False
    return True








class CommandLine:
    def __init__(self):
        self.snakemake="snakemake --rerun-incomplete -k "
        self.runid="".join(random.choices("abcdefghisz",k=3) + random.choices("123456789",k=5))
        self.config=[]
        self.configfile_loaded=False
        self.is_segout_sample=False
        self.is_this_an_segout_run=False
        self.parameters=dict()
        self.log=True
        
    def __str__(self):
        return self.snakemake
    def __repr__(self):
        return self.snakemake
    

    def add_config_argument(self):
        self.snakemake = self.snakemake + " --config " + " ".join(self.config)


    def load_configfile_if_available(self,arguments):
        if self.configfile_loaded is False:
            if "--configfile" in arguments and arguments["--configfile"]:
                self.snakemake = self.snakemake + " --configfile={}".format(arguments["--configfile"])
                configfile=arguments["--configfile"]
            else:
               self.snakemake = self.snakemake + " --configfile={}".format(spatialsnake_path + "/config.yaml")
               configfile=spatialsnake_path + "/config.yaml"
            with open(configfile) as f:
               self.parameters=yaml.load(f,Loader=SafeLoader)
            self.configfile_loaded=True

    def prepare_arguments(self,arguments):
        jobs = arguments.get('--jobs', 4)
        self.snakemake = self.snakemake + " -j {} ".format(jobs)
        self.snakemake = self.snakemake +  " -s {} ".format(f"{spatialsnake_path}/workflow/Snakefile")
        self.load_configfile_if_available(arguments)
        if self.is_this_an_segout_run is False:
            self.config.append("datafolder={}".format(arguments['<INPUT>']))
        if "annotion" in arguments:
            self.config.append("annotion_list={}".format(arguments['<INPUT>']))
        else:
            self.config.append("sample_list={}".format(arguments['<INPUT>']))
        self.config.append(f"spatialsnake_path={spatialsnake_path}/")
        for i,b in arguments.items():
            if i not in ["--jobs","--configfile","--option","--unlock","--remove","--dry","--help","--version","<INPUT>","<command>","--install-packages","<TYPE>","segout","<INTEGRATED_ZARR>"]:
                k=i.lstrip("--")
                if k in ["min_cells", "min_genes", "x1", "x2", "y1", "y2", "workers", "threads"]:
                  try:
                    int(b)
                  except ValueError:
                    print(f"Error: {k} must be an integer")
                    sys.exit(1)
                if self.configfile_loaded is False: #if there is no config file, add all parameters given by the command line or defaults. command line parameters have priority over config file parameters
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
        if self.is_this_an_segout_run is False:
            self.config.append("channel={}".format(arguments['<command>']))
            self.config.append("type={}".format(arguments["<TYPE>"]))
        if self.is_this_an_segout_run:
            self.config.append("channel=integration")
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
        
    
    def write_to_log(self,start):
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        logname = f"spatialsnake_{self.runid}_{timestamp}_runlog.log"
        stop = timeit.default_timer()
        if self.log:
            with open(logname,"w") as f:
                f.write(__logo__ + "\n")
                f.write("Run ID : " + self.runid + "\n")
                f.write("spatialsnake version : " + __version__ + "\n")
                f.write("spatialsnake arguments : " + " ".join(sys.argv) + "\n\n")
                f.write("------------------------------" + "\n")
                f.write("Snakemake arguments : " + str(self.snakemake) + "\n\n")
                f.write("------------------------------" + "\n")
                f.write("Run parameters:\n")
                for key, value in sorted(self.parameters.items()):
                    if value is not None and value != "":
                        f.write(f"  {key.ljust(25)} {value}\n")
                f.write("\n")
                f.write("Total run time: {t:.2f} mins \n".format(t=(stop-start)/60))
                f.write("Useful Information:\n")
                f.write("-" * 20 + "\n")
                f.write("  - For help: spatialsnake --help\n")
                f.write("  - To view config file to learn detail params statement: spatialsnake produce-file\n")
                f.write("  - Output files are stored in the 'results' directory\n")
                f.write("  - Use 【spatialsnake [command] [..] --config-file config.yaml】 to run spatialsnake if you want to customize the parameter settings in config.yaml file. \n")
                f.write("=" * 60 + "\n")
                
                
                
                
                
                
                
                
                
                
                
                
                
                
def run_segout(arguments):
    start = timeit.default_timer()
    snakemake_argument=CommandLine()
    snakemake_argument.is_this_an_segout_run = True
    snakemake_argument.prepare_arguments(arguments)
    subprocess.check_call(str(snakemake_argument),shell=True)
    snakemake_argument.write_to_log(start)

def run_workflow(arguments):
    start = timeit.default_timer()
    snakemake_argument=CommandLine()
    if "segout" in arguments and arguments["segout"]:
        snakemake_argument.is_segout_sample = True
    snakemake_argument.prepare_arguments(arguments)
    subprocess.check_call(str(snakemake_argument),shell=True)
    snakemake_argument.write_to_log(start)


def main():
        cli_arguments = docopt(__doc__, version=__version__)
        if cli_arguments["produce-file"]:
            print("Generating config.yaml file and sample.txt and annotion.txt filter.txt")
            print("You can use this as a template for a spatialsnake run. You may change the settings.")
            shutil.copyfile(spatialsnake_path + "/config.yaml", 'config.yaml')
            print("Generating metadata.csv file...")
            with open("sample.txt","w") as f:
                f.write("sample_id\tpath_to_dir\tbin_size\tgroup\n")
            with open("annotion.txt","w") as f:
                f.write("sample\t0\t1\t2\t3\t4\t5.........please input anno by order of cluster\n")
            with open("filter_list.txt","w") as f:
                f.write("sample_id\tmin_cells\tmin_genes\n")    
            return
        if cli_arguments["install-packages"]:
            subprocess.check_call(spatialsnake_path + "/workflow/scripts/install_packages.R")
            return
        
        if not check_command_line_arguments(cli_arguments):
            print("""Please check your command line arguments. Use "cellsnake --help" for more information""")
            return

        print(cli_arguments['<command>'])
        if cli_arguments['<command>'] == 'single_analysis':
            run_workflow(cli_arguments)
        if cli_arguments['<command>'] == 'compare_analysis':
            run_workflow(cli_arguments)
        if cli_arguments['<command>'] == 'segout':
            run_segout(cli_arguments)
        if cli_arguments['<command>'] == 'spatial_ATAC':
            run_workflow(cli_arguments)
                
                
                
                
                
