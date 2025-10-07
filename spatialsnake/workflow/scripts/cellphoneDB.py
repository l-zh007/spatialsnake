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
parser.add_argument('--threshold', type=float, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--threads', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--pvalue', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--output_name', type=str, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()










output_dir=os.path.dirname(args.output_zarr_path)

display(HTML(db_releases_utils.get_remote_database_versions_html()['db_releases_html_table']))
cpdb_version ='v5.0.0'
cpdb_target_dir = os.path.join(output_dir,'cellphonedb_v500_NatProtocol/', cpdb_version)

if not os.path.isdir(cpdb_target_dir):
  db_utils.download_database(cpdb_target_dir, cpdb_version)


if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
    adata.obs["cell_id"]=adata.obs.index
adata.write(args.output_zarr_path)
df_extract = adata.obs[['cell_id', 'celltype']].copy()         ############spot_id > cell_id updata
txt_path = os.path.join(output_dir, f"{args.sample_id}_cellid_cell_type.txt")
df_extract.to_csv(txt_path, sep="\t", index=False)
print(f"已保存两列数据：{txt_path}")

print(adata.obs)



# if channel=="compare_analysis":
#   if type=="slide_seq":
#     print("secceed")
#   unique_samples = adata.obs['region'].unique()
#   for sample in unique_samples:
#       adata_sample = adata[adata.obs['region'] == sample, :].copy()
#       h5ad_path = os.path.join(output_dir, f"{region}.h5ad")
#       adata_sample.write_h5ad(h5ad_path)
#       df_extract = adata_sample.obs[['cell_id', key]].copy()
#       txt_path = os.path.join(output_dir, f"{region}_cellid_cell_type.txt")
#       df_extract.to_csv(txt_path, sep="\t", index=False)
#       print(f"已保存两列数据：{txt_path}")
#   exit()











cpdb_file_path = os.path.expanduser(os.path.join(cpdb_target_dir,"cellphonedb.zip"))
meta_file_path = os.path.join(output_dir, f"{args.sample_id}_cellid_cell_type.txt")
counts_file_path = args.output_zarr_path
out_path = os.path.join(output_dir,f"cellphonedb_output_{args.output_name}")

os.makedirs(out_path, exist_ok=True)



cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_file_path,
    meta_file_path=meta_file_path,
    counts_file_path=counts_file_path,
    counts_data='hgnc_symbol',
    iterations=500,
    threshold=0.1,
    threads=32,
    pvalue=0.05,
    # active_tfs_file_path=active_tf_path,    # 如果使用转录因子分析
    # microenvs_file_path=microenvs_file_path, # 如果使用微环境定义
    # subsampling=True,                     # 大数据集时可启用子抽样
    output_path=out_path,
    output_suffix=args.output_name)

print("CellPhoneDB 分析已完成！")
print(f"结果保存在: {out_path}")





