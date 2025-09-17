import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "1"
import spatialdata as spd
# import spatialdata_plot as splt
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
from sklearn.cluster import KMeans
############################clustering################################
type="xenium"
RES = 0.5 # clustering resolution
NEIGHBORS = 30  # number of neighbors

MIN_DIST=0.5 #default 0.5
SPREAD=2 #default 1
harmony=True
output_dir = "clustering_results"
os.makedirs(output_dir, exist_ok=True)

# samples = {
#     "Kidney_Cancer":["./data/Kidney_Cancer_data","Kidney_Cancer.zarr"],
#     "Kidney_Normal":["./data/Kidney_Normal_data","Kidney_Normal.zarr"]}
# samples = {"Colon_Cancer_P1":["data/Visium_HD_Human_Colon_Cancer_P1",8,"Colon_Cancer_P1.zarr"],
#             "Colon_Cancer_P2":["data/Visium_HD_Human_Colon_Cancer_P2",8,"Colon_Cancer_P2.zarr"]}
samples = {"Non_Lesional_1":["./data/ST_21_NL","Non_Lesional_1.zarr"],
          "Non_Lesional_2":["./data/ST_22_NL","Non_Lesional_2.zarr"]}
pre_dir="../preprocess/preprocessing_results"

for key, inputs in samples.items():
    if type=="slide_seq":
      print(os.path.join(pre_dir,inputs[-1]))
      adata = sc.read_h5ad(os.path.join(pre_dir,inputs[-1]))
    else:
      concatenated_sdata = spd.read_zarr(os.path.join(pre_dir,inputs[-1]))
      for table in concatenated_sdata.tables.keys():
        print(table)
        adata = concatenated_sdata[table]
    # sc.tl.louvain(adata, resolution=RES,key_added="clusters")
    sc.tl.leiden(adata, flavor="igraph",key_added="clusters", resolution=RES,random_state=0)
    

    # kmeans = KMeans(n_clusters=18, random_state=0)
    # adata.obs['clusters'] = kmeans.fit_predict(adata.obsm['X_pca'][:, :40])
    # adata.obs['clusters'] = adata.obs['clusters'].astype('category')
    sc.pl.umap(adata,color=[
        "total_counts",
        "n_genes_by_counts",
        "clusters"],wspace=0.4)
    plt.savefig(
        os.path.join(output_dir, f"{key}UMAP.png"),
        dpi=300,
        bbox_inches='tight')
    plt.show()
    plt.close()
    # sc.pl.tsne(adata, color=[
    #     "total_counts",
    #     "n_genes_by_counts",
    #     "clusters"], ncols=3, show=True)
    # plt.savefig(
    #     os.path.join(output_dir, "tsene.png"),
    #     dpi=300,
    #     bbox_inches='tight')
    # plt.show()
    # plt.close()
    
    sample_names = adata.obs["region"].unique()
    plt.imshow(pd.crosstab(adata.obs["region"], adata.obs["clusters"]), cmap='hot', interpolation='nearest')
    plt.title("Cell Distribution Across Clusters")
    plt.xlabel("Cluster")
    plt.yticks(range(len(sample_names)), sample_names)
    plt.savefig(
        os.path.join(output_dir, f"{key}Cell_Distribution_Across_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
    plt.show()
    plt.close()
    if type!="slide_seq":
      concatenated_sdata[table]=adata
      concatenated_sdata.write(
        os.path.join(output_dir, f"{inputs[-1]}"),
        overwrite=True)
    else:  
      adata.write(os.path.join(output_dir, f"{inputs[-1]}"))
