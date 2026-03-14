import os
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
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
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir in sys.path:
  sys.path.remove(script_dir)
import squidpy as sq
import cellcharter as cc
import scvi
from lightning.pytorch import seed_everything
from spatialsnake.workflow.function.plot import cluster_proportion
from scipy import sparse as sp
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv

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
parser.add_argument('--channal', type=str, required=True,
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

type = args.type
dir_path = os.path.dirname(args.output_zarr_path)
channel = args.channal
sample_id = args.sample_id

def has_negative(x):
  if x is None:
    return False
  if sp.issparse(x):
    return x.data.size > 0 and x.data.min() < 0
  return np.min(x) < 0

if args.cellcharter_col is None:
  args.cellcharter_col = "spatial_cluster"
if args.max_cluster is None:
  args.max_cluster = 10
if args.significance is None:
  args.significance = 0.05

concatenated_sdata = None
table_key = None
images = []
shapes = []
systems = []

if type == "slide_seq" or args.input_dir.endswith(".h5ad"):
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  table_key = next(iter(concatenated_sdata.tables.keys()))
  image_elements = list(concatenated_sdata.images.keys())
  shape_elements = list(concatenated_sdata.shapes.keys())
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  if len(valid_coord_systems) > 1:
    systems = [s for s in valid_coord_systems if args.image_type in s]
    if len(systems) == 0:
      systems = valid_coord_systems
  else:
    systems = valid_coord_systems
  if len(image_elements) > 1:
    images = [i for i in image_elements if args.image_type in i]
    if len(images) == 0:
      images = image_elements
  else:
    images = image_elements
  if len(shape_elements) > 1 and type != "visium":
    shapes = [s for s in shape_elements if args.shape_type in s]
    if len(shapes) == 0:
      shapes = shape_elements
  else:
    shapes = shape_elements

  if channel == "compare_analysis" and len(systems) > 1:
    adata_list = []
    for s1 in systems:
      data = to_legacy_anndata(concatenated_sdata, include_images=False, coordinate_system=s1)
      if args.sample_col and args.sample_col not in data.obs.columns:
        data.obs[args.sample_col] = s1
      data.obsm['spatial'][:, 1] = np.max(data.obsm['spatial'][:, 1]) - data.obsm['spatial'][:, 1]
      adata_list.append(data)
    adata = ad.concat(adata_list, axis=0, merge='same', pairwise=True, index_unique='_')
  else:
    coord_system = systems[0] if len(systems) > 0 else None
    adata = to_legacy_anndata(concatenated_sdata, include_images=True, coordinate_system=coord_system)

counts_source = None
if adata.raw is not None:
  raw_x = adata.raw.X
  if not has_negative(raw_x):
    counts_source = raw_x
if counts_source is None and "counts" in adata.layers and not has_negative(adata.layers["counts"]):
  counts_source = adata.layers["counts"]
if counts_source is None and not has_negative(adata.X):
  counts_source = adata.X
if counts_source is None:
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  counts_source = adata.X
adata.layers["counts"] = counts_source

if channel == "compare_analysis":
  if args.sample_col and args.sample_col not in adata.obs.columns and "library_id" in adata.obs.columns:
    adata.obs[args.sample_col] = adata.obs["library_id"]
  if args.sample_col and args.sample_col in adata.obs.columns:
    adata.obs[args.sample_col] = pd.Categorical(adata.obs[args.sample_col])
    adata.uns["spatial"] = {s: {} for s in adata.obs[args.sample_col].unique()}

if channel == "compare_analysis" and args.sample_col:
  scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.sample_col)
else:
  scvi.model.SCVI.setup_anndata(adata, layer="counts")

print(adata)
print("@@@@@@@@@")


model = scvi.model.SCVI(adata)
model.train(early_stopping=True, enable_progress_bar=True)
adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

if channel == "compare_analysis" and args.sample_col:
  sq.gr.spatial_neighbors(adata, library_key=args.sample_col, coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)
  cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter', sample_key=args.sample_col)
else:
  sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)
  cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter')

autok = cc.tl.ClusterAutoK(
  n_clusters=(2, args.max_cluster), 
  max_runs=10,
  convergence_tol=0.001)
autok.fit(adata, use_rep='X_cellcharter')
cc.pl.autok_stability(autok)
plt.savefig(
  os.path.join(dir_path, f"{sample_id}_celchar.png"),
  dpi=300,
  bbox_inches='tight')
plt.close()
adata.obs[args.cellcharter_col] = autok.predict(adata, use_rep='X_cellcharter')
if concatenated_sdata is not None and table_key is not None:
  concatenated_sdata[table_key] = adata

if args.celltype_col and args.celltype_col in adata.obs.columns:
  cc.gr.enrichment(adata, group_key=args.cellcharter_col, label_key=args.celltype_col)
  cc.pl.enrichment(adata, group_key=args.cellcharter_col, label_key=args.celltype_col, figsize=(4,4), fontsize=6, dot_scale=1)
  plt.savefig(
    os.path.join(dir_path, f"{sample_id}_enrichment.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

condition_groups = []
if channel == "compare_analysis" and args.condition_col and args.condition_col in adata.obs.columns:
  condition_groups = list(pd.unique(adata.obs[args.condition_col]))
  for cond in condition_groups:
    adata_cond = adata[adata.obs[args.condition_col] == cond].copy()
    cc.gr.nhood_enrichment(
      adata_cond,
      cluster_key=args.cellcharter_col)
    cc.pl.nhood_enrichment(
      adata_cond,
      cluster_key=args.cellcharter_col,
      annotate=True,
      vmin=-1,
      vmax=1,
      figsize=(3,3),
      fontsize=5)
    plt.savefig(
      os.path.join(dir_path, f"{sample_id}_{cond}_enrichment.png"),
      dpi=300,
      bbox_inches='tight')
    plt.close()
    if args.sample_col and args.sample_col in adata_cond.obs.columns:
      cluster_proportion(adata_cond, sample_col=args.sample_col, 
                        cluster_col=args.cellcharter_col,
                        figsize=(12, 8),
                        palette='tab20',
                        sort_samples=None,
                        sort_clusters=None,
                        title="cellcharter_celltype_proportion",
                        save_path=os.path.join(dir_path, f"{sample_id}_Clusters_proportion.png"),
                        dpi=300)
else:
  cc.gr.nhood_enrichment(
    adata,
    cluster_key=args.cellcharter_col)
  cc.pl.nhood_enrichment(
    adata,
    cluster_key=args.cellcharter_col,
    annotate=True,
    vmin=-1,
    vmax=1,
    figsize=(3,3),
    fontsize=5)
  plt.savefig(
    os.path.join(dir_path, f"{sample_id}_nhood_enrichment.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

if channel == "compare_analysis" and len(condition_groups) > 0 and args.sample_col and args.sample_col in adata.obs.columns:
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
    condition_groups=condition_groups,
    annotate=True,
    figsize=(3,3),
    significance=args.significance,
    fontsize=5)
  plt.savefig(
    os.path.join(dir_path, f"{sample_id}_diff_enrichment.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

if concatenated_sdata is not None and len(images) > 0 and len(shapes) > 0:
  if len(images) != len(shapes):
    print("Check the spatial data to make sure that for every image there is a shape")
    print(images, shapes)
    exit()
  for i in range(len(images)):
    sys = systems[0] if len(systems) < 2 else systems[i]
    title = f"{sample_id}_{images[i]}"
    axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
    concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image", coordinate_systems=sys)
    concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i], color=args.cellcharter_col).pl.show(ax=axes[1], coordinate_systems=sys, title=title)
    plt.savefig(
      os.path.join(dir_path, f"{sample_id}_{images[i]}_Clusters.png"),
      dpi=300,
      bbox_inches='tight')
    plt.close()

export_cluster_csv(
  concatenated_sdata if concatenated_sdata is not None else adata,
  args.type,
  dir_path,
  cell_id_col="cell_id",
  info_col=args.cellcharter_col,
  sample_col=args.sample_col if args.sample_col else "sample",
  sample_id=args.sample_id
)

if concatenated_sdata is not None:
  concatenated_sdata.write(args.output_zarr_path)
else:
  if args.output_zarr_path.endswith(".h5ad"):
    adata.write(args.output_zarr_path)
  else:
    adata.write_zarr(args.output_zarr_path)
