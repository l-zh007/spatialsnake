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
import zarr
import matplotlib.pyplot as plt
import seaborn as sns
output_dir = "qc_results"
os.makedirs(output_dir, exist_ok=True)
type="visium_HD"

# samples = {
#     "normal_1": ["data/normal1.h5ad","normal_1.h5ad"],
#     "normal_2": ["data/normal2.h5ad","normal_2.h5ad"]
# }
samples = {"Colon_Cancer_P1":["data/Visium_HD_Human_Colon_Cancer_P1",8,"Colon_Cancer_P1.zarr"],
            "Colon_Cancer_P2":["data/Visium_HD_Human_Colon_Cancer_P2",8,"Colon_Cancer_P2.zarr"]}
for key, inputs in samples.items():
  if type=="slide_seq":
    adata = sc.read_h5ad(inputs[-1])
    adata.obs['region']=key
  else:
    concatenated_sdata = spd.read_zarr(inputs[-1])
    for table in concatenated_sdata.tables.keys():
      print(table)
      adata = concatenated_sdata[table]
  print(adata)
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))


  sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
#### visium wu
#控制探针和控制代码词的百分比可以从 adata.obs 中计算得出。
# cprobes = (
#     concatenated_sdata["segmentation_counts"].obs["control_probe_counts"].sum() / concatenated_sdata["segmentation_counts"].obs["total_counts"].sum() * 100
# )
# cwords = (
#     concatenated_sdata["segmentation_counts"].obs["control_codeword_counts"].sum() / concatenated_sdata["segmentation_counts"].obs["total_counts"].sum() * 100
# )
# print(f"Negative DNA probe count % : {cprobes}")
# print(f"Negative decoding count % : {cwords}")



  fig, axs = plt.subplots(1, 2, figsize=(15, 4))
  axs[0].set_title("Total transcripts per cell")
  sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0])

  axs[1].set_title("Unique transcripts per cell")
  sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1])

#### visium wu
# axs[2].set_title("Area of segmented cells")
# sns.histplot(
#     concatenated_sdata["segmentation_counts"].obs["cell_area"],
#     kde=False,
#     ax=axs[2],
# )
# 
# axs[3].set_title("Nucleus ratio")
# sns.histplot(
#     concatenated_sdata["segmentation_counts"].obs["nucleus_area"] / concatenated_sdata["segmentation_counts"].obs["cell_area"],
#     kde=False,
#     ax=axs[3],
# )

  plt.savefig(
      os.path.join(output_dir, "total.png"),
      dpi=300,
      bbox_inches='tight')
  plt.close()









  sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=['mt'], 
    inplace=True, 
    percent_top=None)




  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts"], 
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total UMI by Sample")
  plt.axhline(y=4, color='r', linestyle='-')
  plt.axhline(y=8, color='r', linestyle='-')
  plt.savefig(
    os.path.join(output_dir, "total_umi_by_sample.png"),
    dpi=300, 
    bbox_inches='tight')
  plt.show()
  plt.close()



  sc.pl.violin(
    adata=adata, 
    keys=["log1p_n_genes_by_counts"], 
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total Genes by Sample")
  plt.savefig(
    os.path.join(output_dir, "total_genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()

  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts_mt"], 
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Mitochondrial Genes by Sample")
  plt.savefig(
    os.path.join(output_dir, "genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()


#################### xenium  #####################

# for table in concatenated_sdata.tables.values():
#         table.obs['spot_id'] = table.obs['spot_id'].astype(str)          ############## spot_id cell_id     
#         table.obs['region'] = table.obs['region'].astype('category')


  if type!="slide_seq":
    concatenated_sdata[table]=adata
    concatenated_sdata.write(
      os.path.join(output_dir, f"{inputs[-1]}"),
      overwrite=True)
  else:  
    adata.write(os.path.join(output_dir, f"{inputs[-1]}"))
