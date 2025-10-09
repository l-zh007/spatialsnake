###################preprocessing###########################################
import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "8"
import sys
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import spatialdata_io
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
import argparse
import seaborn as sns



parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')                   

parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')

parser.add_argument('--min_genes', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--min_cells', type=int, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--variable', type=bool, required=False,
                   help='Path for the output zarr file')                    
parser.add_argument('--filter_dict',type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--seg_filter', type=bool, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--NEIGHBORS', type=int, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--harmony', type=bool, required=True,
                   help='Path for the output zarr file')                
args = parser.parse_args()
type=args.type
dir_path=os.path.dirname(args.output_zarr_path)



if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]



if args.filter_dict:
  print(args.filter_dict,args.sample_id)
  filter_dict = json.loads(args.filter_dict)
  min_cells=filter_dict[args.sample_id][0]
  min_genes=filter_dict[args.sample_id][1]
  
  
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
  adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
  adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
  sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"], 
    percent_top=(10, 20, 50, 150),
    inplace=True, 
    log1p=True)
else:
  min_genes=args.min_genes
  min_cells=args.min_cells
print(min_cells,min_genes)

# if seg_filter==True:
#   adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
#   adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
#   adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
#   sc.pp.calculate_qc_metrics(
#     adata, 
#     qc_vars=["mt", "ribo", "hb"], 
#     percent_top=(10, 20, 50, 150),
#     inplace=True, 
#     log1p=True)


sc.pp.filter_genes(adata,min_cells=min_cells)
sc.pp.filter_cells(adata,min_counts=min_genes)

# sc.pp.filter_cells(adata, min_genes=min_genes)

sc.pl.violin(adata=adata, keys=["log1p_total_counts"], stripplot=False, inner="box",show=False, groupby="region",)
plt.title("Total UMI by Sample")
plt.axhline(y=4, color='r', linestyle='-')
plt.axhline(y=8, color='r', linestyle='-')
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}filtered_Total_UMI.png"),
  dpi=300,
  bbox_inches='tight')
plt.show()
plt.close()

sc.pl.violin(adata=adata, keys=["log1p_n_genes_by_counts"], groupby="region", stripplot=False, inner="box",show=False)
plt.title("Total Genes by Sample")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}filtered_Total_Genes.png"),
    dpi=300,
    bbox_inches='tight')
plt.show()
plt.close()

sc.pl.violin(adata=adata, keys=["log1p_total_counts_mt"], groupby="region", stripplot=False, inner="box",show=False)
plt.title("Mitochondrial Genes by Sample")
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}_Mitochondrial_Genes.png"),
  dpi=300,
  bbox_inches='tight')
plt.show()
plt.close()

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_scatter.png"),
    dpi=300,
    bbox_inches='tight')
plt.show()
plt.close()


sc.pp.normalize_total(adata, target_sum = None)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# print(f"高变异基因数量: {sum(adata.var.highly_variable)}")
sc.pl.highly_variable_genes(adata)
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}highly_variable.png"),
  dpi=300,
  bbox_inches='tight')
plt.show()
plt.close()

if args.variable==True:
  adata = adata[:, adata.var.highly_variable].copy()



sc.tl.pca(adata)


if args.harmony:
    sce.pp.harmony_integrate(adata, key="region", basis="X_pca",max_iter_harmony=20)

sc.pp.neighbors(adata, n_neighbors=args.NEIGHBORS, n_pcs=40)

  
sc.pl.pca_variance_ratio(adata, log=True,n_pcs=50)
plt.title("pca_variance_ratio")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}pca_variance_ratio.png"),
    dpi=300,
    bbox_inches='tight')
plt.show()
plt.close()
  





if type!="slide_seq":
  adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
  adata.obs['region'] = adata.obs['region'].astype('category')
  concatenated_sdata[table]=adata
  # for table in concatenated_sdata.tables.values():
  #     table.obs['spot_id'] = table.obs['cell_id'].astype(str)
  #     table.obs['region'] = table.obs['region'].astype('category')
  concatenated_sdata.write(os.path.join(args.output_zarr_path),overwrite=True)
else:
  adata.write(os.path.join(args.output_zarr_path))

del concatenated_sdata,adata
gc.collect()
