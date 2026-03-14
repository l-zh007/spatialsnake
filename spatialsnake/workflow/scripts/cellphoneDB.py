import pandas as pd
import glob
import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "1"
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils
from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

import argparse
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--counts_data', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--iterations', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--threads', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--pvalue', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--output_name', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--microenvs_file_path', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--active_tf_path', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--threshold', type=float, required=False,
                   help='Path for the output zarr file')   #0.1                
parser.add_argument('--degs_file_path', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--celltype_col', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--niche_col', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--is_singlecell', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--method', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--de_method', type=str, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()
output_dir=os.path.dirname(args.output_zarr_path)
display(HTML(db_releases_utils.get_remote_database_versions_html()['db_releases_html_table']))
cpdb_version ='v5.0.0'
cpdb_target_dir = os.path.join(output_dir,'cellphonedb_v500_NatProtocol/', cpdb_version)

if not os.path.isdir(cpdb_target_dir):
  db_utils.download_database(cpdb_target_dir, cpdb_version)

def parse_bool(value):
  if isinstance(value, bool):
    return value
  if value is None:
    return False
  value_str = str(value).strip().lower()
  if value_str in {"true", "1", "yes", "y", "t"}:
    return True
  if value_str in {"false", "0", "no", "n", "f", ""}:
    return False
  raise ValueError(f"Invalid boolean value: {value}")

run_type = args.type
input_path = args.input_dir
celltype_col = args.celltype_col
niche_col = args.niche_col
is_singlecell = parse_bool(args.is_singlecell)
method = args.method if args.method else "statistical"
de_method = args.de_method if args.de_method else "wilcoxon"
iterations = args.iterations if args.iterations is not None else 500
threads = args.threads if args.threads is not None else 32
pvalue = args.pvalue if args.pvalue is not None else 0.05
input_ext = os.path.splitext(input_path)[1].lower()

if input_ext in [".h5ad", ".h5"] or run_type == "slide_seq":
  adata = sc.read_h5ad(input_path)
else:
  concatenated_sdata = spd.read_zarr(input_path)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
    adata.obs["cell_id"]=adata.obs.index

if celltype_col not in adata.obs.columns:
  raise ValueError(f"celltype_col not found in obs: {celltype_col}")

os.makedirs(output_dir, exist_ok=True)
adata.write(args.output_zarr_path)
df_extract = adata.obs[['cell_id', celltype_col]].copy()
df_extract = df_extract.rename(columns={celltype_col: "cell_type"})
txt_path = os.path.join(output_dir, f"{args.sample_id}_cellid_cell_type.txt")
df_extract.to_csv(txt_path, sep="\t", index=False)

cpdb_file_path = os.path.expanduser(os.path.join(cpdb_target_dir,"cellphonedb.zip"))
meta_file_path = txt_path
counts_file_path = args.output_zarr_path
out_path = os.path.join(output_dir,f"cellphonedb_output")
os.makedirs(out_path, exist_ok=True)

microenvs_file_path = None
if args.microenvs_file_path and os.path.isfile(args.microenvs_file_path):
  microenvs_file_path = args.microenvs_file_path
if not is_singlecell and microenvs_file_path is None:
  if not niche_col:
    raise ValueError("niche_col is required for spatial data")
  if niche_col not in adata.obs.columns:
    raise ValueError(f"niche_col not found in obs: {niche_col}")
  microenvs_df = adata.obs[[celltype_col, niche_col]].copy()
  microenvs_df = microenvs_df.dropna()
  microenvs_df[celltype_col] = microenvs_df[celltype_col].astype(str)
  microenvs_df[niche_col] = microenvs_df[niche_col].astype(str)
  microenvs_df = microenvs_df.rename(columns={celltype_col: "cell_type", niche_col: "microenvironment"})
  microenvs_file_path = os.path.join(output_dir, f"{args.sample_id}_microenvs.txt")
  microenvs_df.to_csv(microenvs_file_path, sep="\t", index=False)

active_tf_path = args.active_tf_path if args.active_tf_path and os.path.isfile(args.active_tf_path) else None

from cellphonedb.src.core.methods import cpdb_degs_analysis_method
if method == "degs":
  if not args.degs_file_path or not os.path.isfile(args.degs_file_path):
    raise ValueError("degs_file_path is required when method is degs")
  cpdb_results = cpdb_degs_analysis_method.call(
         cpdb_file_path = cpdb_file_path,
         meta_file_path = meta_file_path,
         counts_file_path = counts_file_path,
         degs_file_path = args.degs_file_path,
         counts_data = args.counts_data,
         threshold = args.threshold,
         output_path = out_path,
         output_suffix=args.output_name)
else:
  cpdb_kwargs = {
    "cpdb_file_path": cpdb_file_path,
    "meta_file_path": meta_file_path,
    "counts_file_path": counts_file_path,
    "counts_data": args.counts_data,
    "iterations": iterations,
    "threshold": args.threshold,
    "threads": threads,
    "pvalue": pvalue,
    "score_interactions": True,
    "output_path": out_path,
    "output_suffix": args.output_name,
  }
  if active_tf_path:
    cpdb_kwargs["active_tfs_file_path"] = active_tf_path
  if microenvs_file_path:
    cpdb_kwargs["microenvs_file_path"] = microenvs_file_path
  cpdb_results = cpdb_statistical_analysis_method.call(**cpdb_kwargs)











