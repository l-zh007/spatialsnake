import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as spd
from spatialdata.models import PointsModel, TableModel
from spatialdata.transformations import get_transformation
import spatialdata_plot
from spatialsnake.workflow.function.plot import cluster_proportion, plot_spatial_clusters
from spatialsnake.workflow.function.stereoseq_selection import (
  extract_table_spatial_metadata,
  resolve_stereoseq_table_key,
  resolve_visual_target,
)
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}
logger = setup_logger("manual_annotation")

def parse_bool(value):
  if isinstance(value, bool):
    return value
  if value is None:
    return False
  value_str = str(value).strip().lower()
  if value_str in ["true", "1", "yes", "y", "t"]:
    return True
  if value_str in ["false", "0", "no", "n", "f", "none", "null", ""]:
    return False
  raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

def select_elements(items, keyword):
    if args.type not in ["visium","visium_HD","visium_segment"]:
      return items
    if not keyword:
        return list(items)
    selected = [item for item in items if keyword in item]
    return selected if selected else list(items)

def pick_table_key(sdata, run_type=None, input_spec=None):
  table_keys = list(getattr(sdata, "tables", {}).keys())
  if not table_keys:
    raise RuntimeError("No table found in SpatialData input.")
  if run_type in STEREOSEQ_TYPES:
    return resolve_stereoseq_table_key(input_spec, table_keys)
  return "table" if "table" in table_keys else table_keys[0]


def parse_annotation_mapping(raw_anno_data, sample_id):
  if not raw_anno_data:
    return {}
  mapping_data = json.loads(raw_anno_data)
  if sample_id in mapping_data and isinstance(mapping_data[sample_id], dict):
    return {str(k): str(v) for k, v in mapping_data[sample_id].items()}
  if "concatenated_sdata" in mapping_data and isinstance(mapping_data["concatenated_sdata"], dict):
    return {str(k): str(v) for k, v in mapping_data["concatenated_sdata"].items()}
  return {}


def annotate_adata(adata, cell_annotation):
  cluster_values = adata.obs["clusters"].astype(str)
  mapped = cluster_values.map(cell_annotation)
  mapped = mapped.fillna("Cluster_" + cluster_values)
  adata.obs["celltype"] = mapped.astype("category")
  return adata


def ensure_obs_columns(adata, sample_id):
  if "cell_id" not in adata.obs.columns:
    adata.obs["cell_id"] = adata.obs_names.astype(str)
  if "region" not in adata.obs.columns:
    adata.obs["region"] = str(sample_id)
  if "group" not in adata.obs.columns:
    adata.obs["group"] = adata.obs["region"].astype(str)
  adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
  adata.obs["region"] = adata.obs["region"].astype(str).astype("category")
  adata.obs["group"] = adata.obs["group"].astype(str)
  return adata


def sanitize_obs_index_name(adata, instance_key):
  if instance_key and instance_key in adata.obs and adata.obs_names.name == instance_key:
    safe_name = "__obs_index__"
    while safe_name in adata.obs:
      safe_name = f"_{safe_name}"
    adata.obs_names.name = safe_name
  return adata


def save_cell_cluster_csv(adata, output_path):
  output_df = pd.DataFrame(index=adata.obs_names.astype(str))
  output_df.index.name = "cell_id"
  output_df["clusters"] = adata.obs["clusters"].astype(str).to_numpy()
  output_df["celltype"] = adata.obs["celltype"].astype(str).to_numpy()
  if "region" in adata.obs.columns:
    output_df["region"] = adata.obs["region"].astype(str).to_numpy()
  if "group" in adata.obs.columns:
    output_df["group"] = adata.obs["group"].astype(str).to_numpy()
  output_df.to_csv(output_path)


def crop0(x, crs, x1, x2, y1, y2):
  return spd.bounding_box_query(
    x,
    min_coordinate=[x1, y1],
    max_coordinate=[x2, y2],
    axes=("x", "y"),
    target_coordinate_system=crs,
  )


parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True, help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True, help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True, help='Path for the output zarr file')
parser.add_argument('--coord', type=int, nargs=4, required=False, help='Path for the output zarr file')
parser.add_argument('--anno_data', type=str, required=False, help='Path for the output zarr file')
parser.add_argument('--image_slice', type=parse_bool, required=False, help='Path for the output zarr file')
parser.add_argument('--vis_mode', type=str, required=False, help='Spatial visualization mode: auto, point, or shape')
parser.add_argument('--point_size', type=float, required=False, default=2.5, help='Point size for point-based spatial visualization')
parser.add_argument('--input_spec', type=str, required=False, help='Stereo-seq input mode passed from sample.txt')
args = parser.parse_args()

dir_path = os.path.dirname(args.output_zarr_path)
os.makedirs(dir_path, exist_ok=True)
cell_annotation = parse_annotation_mapping(args.anno_data, args.sample_id)
run_type = args.type
log_step(logger, 1, 5, "loading clustered data")

concatenated_sdata = spd.read_zarr(args.input_dir)
table = pick_table_key(concatenated_sdata, run_type=run_type, input_spec=args.input_spec)
adata = concatenated_sdata[table]
adata = ensure_obs_columns(adata, args.sample_id)
adata = annotate_adata(adata, cell_annotation)
table_meta = extract_table_spatial_metadata(adata)
region_name = table_meta["region_name"]
instance_key = table_meta["instance_key"]
if isinstance(region_name, (list, tuple)):
  region_name = region_name[0] if region_name else None
adata.obs_names = adata.obs_names.astype(str)
adata = sanitize_obs_index_name(adata, instance_key)
concatenated_sdata[table] = adata
#adata.obs[instance_key] = adata.obs_names
log_step(logger, 2, 5, "rendering cell type spatial plots")
plot_spatial_clusters(
  concatenated_sdata, table, "celltype", dir_path, args.sample_id,
  vis_mode=args.vis_mode, point_size=args.point_size,
  coord=args.coord if args.image_slice else None,
)

log_step(logger, 3, 5, "saving cell type proportion plot")
cluster_proportion(
  adata,
  sample_col='region',
  cluster_col='celltype',
  figsize=(12, 8),
  palette='tab20',
  sort_samples=None,
  sort_clusters=None,
  title="celltype_proportion",
  save_path=os.path.join(dir_path, "celltype_proportion.png"),
  dpi=300
)

log_step(logger, 4, 5, "saving UMAP and summary plots")
if "X_umap" in adata.obsm:
  sc.pl.umap(adata, color="celltype", wspace=0.4, show=False)
  plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_UMAP.png"),
    dpi=300,
    bbox_inches='tight'
  )
  plt.close()

if "n_genes_by_counts" in adata.obs.columns:
  sc.pl.violin(adata, ["n_genes_by_counts"], groupby="celltype", rotation=90, show=False)
  plt.ylabel("Detected genes")
  plt.title("Detected genes by annotation")
  plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_annotation_summary.png"),
    dpi=300,
    bbox_inches='tight'
  )
  plt.close()
else:
  counts = adata.obs["celltype"].astype(str).value_counts().sort_values()
  plt.figure(figsize=(10, max(6, len(counts) * 0.3)))
  counts.plot(kind="barh")
  plt.title("Observation counts by annotation")
  plt.xlabel("Observation count")
  plt.tight_layout()
  plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_annotation_summary.png"),
    dpi=300,
    bbox_inches='tight'
  )
  plt.close()

save_cell_cluster_csv(adata, os.path.join(dir_path, f"{args.sample_id}_cell_clusters.csv"))

log_step(logger, 5, 5, f"saving annotated data to {args.output_zarr_path}")
adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
adata.obs[instance_key] = adata.obs[instance_key].astype(str)
adata.obs['region'] = adata.obs['region'].astype('category')
concatenated_sdata[table] = adata
concatenated_sdata.write(os.path.join(args.output_zarr_path), overwrite=True)
logger.info("Manual annotation module completed")
