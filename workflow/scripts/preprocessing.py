###################preprocessing###########################################
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
MIN_DIST=0.5 #default 0.5
SPREAD=2 #default 1
type="xenium"
tsene=False
samples = {"Non_Lesional_1":["./data/ST_21_NL","Non_Lesional_1.zarr"],
          "Non_Lesional_2":["./data/ST_22_NL","Non_Lesional_2.zarr"]}

# samples = {
#     "normal_1": ["data/normal1.h5ad","normal_1.h5ad"],
#     "normal_2": ["data/normal2.h5ad","normal_2.h5ad"]
# }
# samples = {
#     "Kidney_Cancer":["./data/Kidney_Cancer_data","Kidney_Cancer.zarr"],
#     "Kidney_Normal":["./data/Kidney_Normal_data","Kidney_Normal.zarr"]}
# samples = {"Colon_Cancer_P1":["data/Visium_HD_Human_Colon_Cancer_P1",8,"Colon_Cancer_P1.zarr"],
#             "Colon_Cancer_P2":["data/Visium_HD_Human_Colon_Cancer_P2",8,"Colon_Cancer_P2.zarr"]}
output_dir = "preprocessing_results"
os.makedirs(output_dir, exist_ok=True)

pre_dir="../intergrate/qc_results"

for key, inputs in samples.items():
  if type=="slide_seq":
    print(os.path.join(pre_dir,inputs[-1]))
    adata = sc.read_h5ad(os.path.join(pre_dir,inputs[-1]))
  else:
    concatenated_sdata = spd.read_zarr(os.path.join(pre_dir,inputs[-1]))
    for table in concatenated_sdata.tables.keys():
      print(table)
      adata = concatenated_sdata[table]

  sc.pp.filter_cells(adata, min_counts=10)
  sc.pp.filter_genes(adata, min_cells=5)


  sc.pl.violin(adata=adata, keys=["log1p_total_counts"], stripplot=False, inner="box",show=False, groupby="region",)
  plt.title("Total UMI by Sample")
  plt.axhline(y=4, color='r', linestyle='-')
  plt.axhline(y=8, color='r', linestyle='-')
  plt.savefig(
    os.path.join(output_dir, f"{key}filtered_Total_UMI.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()

  sc.pl.violin(adata=adata, keys=["log1p_n_genes_by_counts"], groupby="region", stripplot=False, inner="box",show=False)
  plt.title("Total Genes by Sample")
  plt.savefig(
    os.path.join(output_dir, f"{key}filtered_Total_Genes.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()

  sc.pl.violin(adata=adata, keys=["log1p_total_counts_mt"], groupby="region", stripplot=False, inner="box",show=False)
  plt.title("Mitochondrial Genes by Sample")
  plt.savefig(
    os.path.join(output_dir, f"{key}_Mitochondrial_Genes.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()


  sc.pp.normalize_total(adata, target_sum = None)
  sc.pp.log1p(adata)
  sc.tl.pca(adata)
  sc.pp.scale(adata, max_value=10)
  sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
  sc.tl.umap(adata,min_dist=MIN_DIST, spread=SPREAD, random_state=0)
  
  if tsene:
    sc.tl.tsne(adata, n_pcs=50, perplexity=30)
  

  sc.pl.pca_variance_ratio(adata, log=True,n_pcs=50)
  plt.title("pca_variance_ratio")
  plt.savefig(
    os.path.join(output_dir, f"{key}pca_variance_ratio.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()



# for table in concatenated_sdata.tables.values():
#         table.obs['spot_id'] = table.obs['spot_id'].astype(str)
#         table.obs['region'] = table.obs['region'].astype('category')
  if type!="slide_seq":
    concatenated_sdata[table]=adata
    concatenated_sdata.write(
      os.path.join(output_dir, f"{inputs[-1]}"),
      overwrite=True)
  else:  
    adata.write(os.path.join(output_dir, f"{inputs[-1]}"))

