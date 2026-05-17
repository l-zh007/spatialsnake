import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "16"
import warnings
warnings.filterwarnings("ignore")
import spatialdata as spd
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import json
import gc
import geopandas as gpd
import matplotlib.pyplot as plt
from PIL import Image
from spatialdata.models import TableModel
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
from sklearn.cluster import KMeans
import argparse
from spatialsnake.workflow.function.stereoseq_selection import extract_table_spatial_metadata, resolve_stereoseq_table_key
############################clustering################################
STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}


def pick_table_key(sdata, run_type=None, input_spec=None):
  table_keys = list(getattr(sdata, "tables", {}).keys())
  if not table_keys:
    raise RuntimeError("No table found in SpatialData input.")
  if run_type in STEREOSEQ_TYPES:
    return resolve_stereoseq_table_key(input_spec, table_keys)
  return "table" if "table" in table_keys else table_keys[0]

parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--tsene', type=str, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--MIN_DIST', type=float, required=False,
                   help='Path for the output zarr file')                    
parser.add_argument('--SPREAD', type=float, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--RES', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--cluster_algorithm', type=str, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--n_clusters', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sketch', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--pcs', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--NEIGHBORS', type=int, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--input_spec', type=str, required=False,
                   help='Stereo-seq input mode passed from sample.txt')
args = parser.parse_args()
type=args.type

sc.settings.n_jobs = 8

dir_path=os.path.dirname(args.output_zarr_path)

if type=="slide_seq":
  adata_for_sketch = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  table = pick_table_key(concatenated_sdata, run_type=type, input_spec=args.input_spec)
  adata_for_sketch = concatenated_sdata[table]



pca_selection = args.pcs
print("手动主成分选择")
  
sc.pp.neighbors(adata_for_sketch, n_neighbors=args.NEIGHBORS, use_rep="X_pca", metric="euclidean", n_pcs=pca_selection)


if args.tsene=="True":
  print("running tsene !!!!!!!!!!!!!")
  sc.tl.tsne(adata_for_sketch, n_pcs=pca_selection, perplexity=30)
      
if args.cluster_algorithm=="leiden":    
  sc.tl.leiden(adata_for_sketch, flavor="igraph",key_added="clusters", resolution=args.RES,random_state=0)
elif args.cluster_algorithm=="louvain":   
  sc.tl.louvain(adata_for_sketch, resolution=args.RES,key_added="clusters")
elif args.cluster_algorithm=="kmeans":
  kmeans = KMeans(n_clusters=args.n_clusters, random_state=0)
  adata_for_sketch.obs['clusters'] = kmeans.fit_predict(adata_for_sketch.obsm['X_pca'][:, :40])
  adata_for_sketch.obs['clusters'] = adata_for_sketch.obs['clusters'].astype('category')


sc.tl.umap(adata_for_sketch,min_dist=args.MIN_DIST, spread=args.SPREAD)


sketch = args.sketch=="True"
if sketch:
  print("using sketch !!!!!!!!!!!!!")
  sketch_adata =sc.read_h5ad(os.path.join(os.path.dirname(args.input_dir),"sketch.h5ad"))
  print(sketch_adata)
  chunk_size = 100000
  num_chunks = sketch_adata.n_obs // chunk_size + (1 if sketch_adata.n_obs % chunk_size != 0 else 0)
  adata_query_ingested_list = []
  for i in range(num_chunks):
      start = i * chunk_size
      end = min((i + 1) * chunk_size, sketch_adata.n_obs)
      adata_query_chunk = sketch_adata[start:end].to_memory()
      sc.tl.ingest(adata_query_chunk, adata_for_sketch, obs="clusters")
      adata_query_ingested_list.append(adata_query_chunk)
  adata_query_ingested = sc.concat(adata_query_ingested_list)
  adata_query_ingested.var=adata_query_chunk.var
  adata_query_ingested.varm=adata_for_sketch.varm
  adata_query_ingested.uns=adata_for_sketch.uns
  if "region_colors" in adata_query_ingested.uns.keys():
    adata_query_ingested.uns['region_colors']=adata_query_ingested.uns['region_colors']
  else:
    adata_query_ingested.uns['clusters_colors']=adata_query_ingested.uns['clusters_colors']
else:
  adata_query_ingested=adata_for_sketch
del adata_for_sketch
sc.pl.umap(adata_query_ingested,color=[
    "total_counts",
    "n_genes_by_counts",
    "clusters"],wspace=0.4)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}UMAP.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()



if args.tsene=="True":
  print("running tsene !!!!!!!!!!!!!")
  sc.pl.tsne(adata_query_ingested, color=[
      "total_counts",
      "n_genes_by_counts",
      "clusters"], ncols=3, show=True)
  plt.savefig(
      os.path.join(dir_path,f"{args.sample_id}tsene.png"),
      dpi=300,
      bbox_inches='tight')
  plt.show()
  plt.close()


sample_names = adata_query_ingested.obs["region"].unique()
plt.imshow(pd.crosstab(adata_query_ingested.obs["region"], adata_query_ingested.obs["clusters"]), cmap='hot', interpolation='nearest')
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
    table_meta = extract_table_spatial_metadata(adata_query_ingested)
    region_name = table_meta["region_name"]
    region_key = table_meta["region_key"]
    instance_key = table_meta["instance_key"]
    adata_query_ingested.obs_names = adata_query_ingested.obs_names.astype(str)
    adata_query_ingested.obs[instance_key] = adata_query_ingested.obs[instance_key].astype(str)
    concatenated_sdata[table]=adata_query_ingested
    #if region_name is not None:
      #adata_query_ingested.obs[region_key] = pd.Categorical([region_name] * adata_query_ingested.n_obs)
    #adata_query_ingested.uns.pop("spatialdata_attrs", None)
    #concatenated_sdata[table] = TableModel.parse(
      #adata_query_ingested,
      #region=region_name,
      #region_key=region_key,
      #instance_key=instance_key,
    #)
    concatenated_sdata.write(
      os.path.join(args.output_zarr_path),
      overwrite=True)
else:  
    adata_query_ingested.write(os.path.join(args.output_zarr_path))
