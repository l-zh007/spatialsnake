import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import torch
torch.set_float32_matmul_precision('high')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import gc
import matplotlib.pyplot as plt
import argparse
from spatialdata.datasets import blobs_annotating_element
from spatialdata.transformations import Affine, set_transformation
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
import warnings
warnings.filterwarnings("ignore")
import squidpy as sq
import cellcharter as cc
import scvi
from lightning.pytorch import seed_everything
from spatialsnake.workflow.function.plot import cluster_proportion

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
parser.add_argument('--image_type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--channal', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--shape_type', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--significance', type=float,required=False,
                   help='Path for the output zarr file') #0.05
parser.add_argument('--max_cluster', type=int,required=False,
                   help='Path for the output zarr file') ##10 
                   
parser.add_argument('--condition_col', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_col', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--celltype_col', type=str,required=False,
                   help='Path for the output zarr file')              
parser.add_argument('--cellcharter_col', type=str,required=False,
                   help='Path for the output zarr file')   
                  
args = parser.parse_args()

type=args.type
dir_path=os.path.dirname(args.output_zarr_path)


if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]

sdata.layers["counts"] = sdata.X.copy()


sample_cnt=args.sample_cnt
image_elements = list(concatenated_sdata.images.keys())
shape_elements = list(concatenated_sdata.shapes.keys())
valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
image_counts=len(image_elements)//sample_cnt
shape_count=len(shape_elements)//sample_cnt
shapes=[]
images=[]
systems=[]
if len(valid_coord_systems)>1:
  for i in range(len(image_elements)):
    if args.image_type in image_elements[i]:
      images.append(image_elements[i])
else:
  images=image_elements
if len(valid_coord_systems)>1:
  for i in range(len(valid_coord_systems)):
      if args.image_type in valid_coord_systems[i]:
        systems.append(valid_coord_systems[i])
else:
  systems=valid_coord_systems
  
print(shape_count,shape_elements)
print(type)
if shape_count>1 and type!="visium":
  for i in range(len(shape_elements)):
    if args.shape_type in shape_elements[i]:
      shapes.append(shape_elements[i])
else:
  shapes=shape_elements

if channel=="compare_analysis":
  for i in range(len(systems)):
    s1=systems[i]
    data = to_legacy_anndata(concatenated_sdata, include_images=False, coordinate_system=s1)
    data.obsm['spatial'][:, 1] = np.max(data.obsm['spatial'][:, 1]) - data.obsm['spatial'][:, 1]
    adata_list.append(data)
  adata = ad.concat(adata_list, axis=0, merge='same', pairwise=True, index_unique='_')
  adata.uns['spatial'] = {s: {} for s in adata.obs[args.sample_col].unique()}
  adata.obs[args.sample_col] = pd.Categorical(adata.obs[args.sample_col])
else:
  s1=systems[0]
  adata = to_legacy_anndata(concatenated_sdata, include_images=False, coordinate_system=s1)




scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key=args.sample_col)
model = scvi.model.SCVI(adata)
model.train(early_stopping=True, enable_progress_bar=True)
adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
sq.gr.spatial_neighbors(adata, library_key=args.sample_col, coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)
cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter', sample_key=args.sample_col)
autok = cc.tl.ClusterAutoK(
    n_clusters=(2,args.max_cluster), 
    max_runs=10,
    convergence_tol=0.001)
autok.fit(adata, use_rep='X_cellcharter')
cc.pl.autok_stability(autok)
plt.savefig(
        os.path.join(dir_path,f"{sample_id}_celchar.png"),
        dpi=300,
        bbox_inches='tight')
plt.close()
sdata.obs[args.cellcharter_col] = autok.predict(adata, use_rep='X_cellcharter')
concatenated_sdata[table] = sdata



condition=adata.obs[args.condition_col].unique() if args.channal=="compare_analysis" else ""
for i in condition:
  if args.channal=="compare_analysis":
    adata_balbc = adata[adata.obs[args.condition_col] == i]
  cc.gr.nhood_enrichment(
    adata_balbc,
    cluster_key=args.cellcharter_col)
  print(adata_balbc)
  cc.pl.nhood_enrichment(
    adata_balbc,
    cluster_key=args.cellcharter_col,
    annotate=True,
    vmin=-1,
    vmax=1,
    figsize=(3,3),
    fontsize=5)
  plt.savefig(
        os.path.join(dir_path,f"{sample_id}_{i}_enrichment.png"),
        dpi=300,
        bbox_inches='tight')
  plt.close()
  cluster_proportion(adata_balbc, sample_col=args.cellcharter_col, 
                      cluster_col=args.celltype_col,
                      figsize=(12, 8),
                      palette='tab20',
                      sort_samples=None,
                      sort_clusters=None,
                      title="cellcharter_celltype_proportion",
                      save_path=os.path.join(dir_path, f"{sample_id}_Clusters_proportion.png"),
                      dpi=300)
                      
                      
                      
                      

if args.channal=="compare_analysis":
  cc.gr.diff_nhood_enrichment(
    adata,
    cluster_key=args.cellcharter_col,
    condition_key=args.condition_col,
    library_key=args.sample_col,
    pvalues=True,
    n_jobs=15,
    n_perms=100)
  cc.pl.diff_nhood_enrichment(
    adata,
    cluster_key=args.cellcharter_col,
    condition_key=args.condition_col,
    condition_groups=condition,
    annotate=True,
    figsize=(3,3),
    significance=args.significance,
    fontsize=5)
  plt.savefig(
        os.path.join(dir_path, f"{sample_id}_diff_enrichment.png"),
        dpi=300,
        bbox_inches='tight')
  plt.close()


if len(images) != len(shapes):
    print("Check the spatial data to make sure that for every image there is a shape")
    print(images,shapes)
    exit()
else:
  for i in range(len(images)):
        sys=systems[0] if len(valid_coord_systems)<2 else systems[i]
        print("Plotting: "+ images[i])
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image",coordinate_systems=sys)
        concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i],color=args.cellcharter_col).pl.show(ax=axes[1],coordinate_systems=sys, title=title)
        plt.savefig(
        os.path.join(dir_path, f"{sample_id}_{images[i]}_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
        plt.close()

concatenated_sdata.write(args.output_zarr_path)
