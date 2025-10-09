import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "8"
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
import argparse
############################clustering################################
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')


parser.add_argument('--tsene', type=bool, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--MIN_DIST', type=float, required=False,
                   help='Path for the output zarr file')                    
parser.add_argument('--SPREAD', type=float, required=False,
                   help='Path for the output zarr file') 

parser.add_argument('--RES', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--harmony', type=bool, required=False,
                   help='Path for the output zarr file') 
                   
parser.add_argument('--cluster_algorithm', type=str, required=False,
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
        
        
        
        
sc.tl.umap(adata,min_dist=args.MIN_DIST, spread=args.SPREAD, random_state=0)
if args.tsene==True:
  sc.tl.tsne(adata, n_pcs=50, perplexity=30)
      
if args.cluster_algorithm=="leiden":    
  sc.tl.leiden(adata, flavor="igraph",key_added="clusters", resolution=args.RES,random_state=0)
elif args.cluster_algorithm=="louvain":   
  sc.tl.louvain(adata, resolution=RES,key_added="clusters")
elif args.cluster_algorithm=="kmeans":
  kmeans = KMeans(n_clusters=18, random_state=0)
  adata.obs['clusters'] = kmeans.fit_predict(adata.obsm['X_pca'][:, :40])
  adata.obs['clusters'] = adata.obs['clusters'].astype('category')



sc.pl.umap(adata,color=[
    "total_counts",
    "n_genes_by_counts",
    "clusters"],wspace=0.4)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}UMAP.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()
    

if args.tsene==True:
  sc.pl.tsne(adata, color=[
      "total_counts",
      "n_genes_by_counts",
      "clusters"], ncols=3, show=True)
  plt.savefig(
      os.path.join(dir_path,f"{args.sample_id}tsene.png"),
      dpi=300,
      bbox_inches='tight')
  plt.show()
  plt.close()



    
sample_names = adata.obs["region"].unique()
plt.imshow(pd.crosstab(adata.obs["region"], adata.obs["clusters"]), cmap='hot', interpolation='nearest')
plt.title("Cell Distribution Across Clusters")
plt.xlabel("Cluster")
plt.yticks(range(len(sample_names)), sample_names)
plt.savefig(
      os.path.join(dir_path, f"{args.sample_id}Cell_Distribution_Across_Clusters.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()




if type!="slide_seq":
    # if type=="visium":
    #   adata.obs['sopt_id'] = adata.obs['spot_id'].astype(str)
    #   adata.obs['region'] = adata.obs['region'].astype('category')
    # else:
    adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
    adata.obs['region'] = adata.obs['region'].astype('category')
    concatenated_sdata[table]=adata
    concatenated_sdata.write(
      os.path.join(args.output_zarr_path),
      overwrite=True)
else:  
    adata.write(os.path.join(args.output_zarr_path))
