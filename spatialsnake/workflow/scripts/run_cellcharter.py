import os
import sys
import torch
if hasattr(torch, "set_float32_matmul_precision"):
  torch.set_float32_matmul_precision("high")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import spatialdata as spd
import scanpy as sc
import matplotlib.pyplot as plt
import argparse
from copy import deepcopy
from spatialsnake.workflow.function.legacy_anndata import to_legacy_anndata
import warnings
import squidpy as sq
import cellcharter as cc
import scvi
from scipy import sparse as sp
from spatialdata.models import PointsModel
from spatialdata.transformations import Identity
from spatialsnake.workflow.function.plot import cluster_proportion
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.stereoseq_selection import resolve_visual_target
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step
import argparse
warnings.filterwarnings("ignore")
logger = setup_logger("cellcharter")
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
parser.add_argument('--channel', type=str, required=True,
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
parser.add_argument('--threads', type=int, default=8,
                   help='Threads allocated to CellCharter')
                  
args = parser.parse_args()
if args.threads < 1:
  parser.error("--threads must be >= 1")
torch.set_num_threads(args.threads)

type = args.type
dir_path = os.path.dirname(args.output_zarr_path)
channel = args.channel
sample_id = args.sample_id

RANDOM_SEED = 0
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
  torch.cuda.manual_seed_all(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED


def _as_region_list(region):
  if region is None:
    return []
  if isinstance(region, (list, tuple, np.ndarray, pd.Index)):
    return [str(value) for value in region]
  return [str(region)]


def _matrix_is_nonnegative_integer(x, atol=1e-6, chunk_size=1_000_000):
  """Return whether a matrix is suitable as an scVI count matrix."""
  if x is None:
    return False
  values = x.data if sp.issparse(x) else np.asarray(x).reshape(-1)
  if values.size == 0:
    return True
  for start in range(0, values.size, chunk_size):
    chunk = np.asarray(values[start:start + chunk_size])
    if not np.isfinite(chunk).all():
      return False
    if np.any(chunk < 0) or np.any(np.abs(chunk - np.rint(chunk)) > atol):
      return False
  return True


def _copy_matrix(x):
  return x.copy() if hasattr(x, "copy") else np.array(x, copy=True)


def _select_counts(adata):
  """Select validated raw counts without accepting log-normalized matrices."""
  if "counts" in adata.layers:
    if _matrix_is_nonnegative_integer(adata.layers["counts"]):
      return _copy_matrix(adata.layers["counts"]), 'layers["counts"]'
    logger.warning('Ignoring layers["counts"] because it is not a non-negative integer matrix')

  if adata.raw is not None:
    raw_var_names = pd.Index(adata.raw.var_names.astype(str))
    var_names = pd.Index(adata.var_names.astype(str))
    if var_names.isin(raw_var_names).all():
      raw_x = adata.raw[:, var_names].X
      if _matrix_is_nonnegative_integer(raw_x):
        return _copy_matrix(raw_x), "raw.X (aligned to adata.var_names)"
      logger.warning("Ignoring raw.X because it is not a non-negative integer matrix")
    else:
      missing = var_names[~var_names.isin(raw_var_names)]
      logger.warning(f"Ignoring raw.X because {len(missing)} adata genes are absent from raw.var_names")

  if _matrix_is_nonnegative_integer(adata.X):
    return _copy_matrix(adata.X), "X"

  raise ValueError(
    'No valid raw counts were found. Expected non-negative, approximately integer values in '
    'layers["counts"], an aligned raw.X, or X; log-normalized values cannot be used as scVI counts.'
  )


def _infer_xy_columns(points_df):
  for x_col, y_col in (("x", "y"), ("center_x", "center_y"),
                       ("coord_x", "coord_y"), ("imagecol", "imagerow")):
    if x_col in points_df.columns and y_col in points_df.columns:
      return x_col, y_col
  raise ValueError("Points element does not contain a supported x/y coordinate pair")


def _validate_spatial_annotations(sdata, table):
  """Validate table-to-element associations, including shape-less Points data."""
  attrs = table.uns.get("spatialdata_attrs", {})
  regions = _as_region_list(attrs.get("region"))
  region_key = attrs.get("region_key")
  instance_key = attrs.get("instance_key")

  if not table.obs_names.is_unique:
    raise ValueError("SpatialData table obs IDs must be unique before CellCharter analysis")
  if not regions:
    raise ValueError("The SpatialData table has no region in uns['spatialdata_attrs']")
  if not instance_key:
    raise ValueError("The SpatialData table has no instance_key in uns['spatialdata_attrs']")

  for region in regions:
    if region not in sdata.shapes and region not in sdata.labels and region not in sdata.points:
      raise ValueError(
        f"Annotated region '{region}' has no matching Shapes, Labels, or Points element"
      )
    if region not in sdata.points:
      continue

    if region_key and region_key in table.obs.columns:
      region_mask = table.obs[region_key].astype(str) == region
      table_obs = table.obs.loc[region_mask].copy()
    elif len(regions) == 1:
      table_obs = table.obs.copy()
    else:
      raise ValueError(
        f"Cannot select table rows for Points region '{region}': region_key '{region_key}' is missing"
      )

    points = sdata.points[region]
    points_df = points.compute() if hasattr(points, "compute") else points.copy()
    points_df = points_df.copy()
    point_index = pd.Index(points_df.index.astype(str))
    if point_index.has_duplicates:
      raise ValueError(f"Points element '{region}' contains duplicate IDs")
    if len(points_df) != len(table_obs):
      raise ValueError(
        f"Points/table row count mismatch for '{region}': {len(points_df)} Points and "
        f"{len(table_obs)} table rows"
      )

    candidate_ids = []
    if instance_key in table_obs.columns:
      table_instance_ids = pd.Index(table_obs[instance_key].astype(str))
      if table_instance_ids.has_duplicates:
        raise ValueError(f"Table instance IDs for Points region '{region}' are not unique")
      candidate_ids.append(("Points index", point_index, table_instance_ids))
      if instance_key in points_df.columns:
        point_instance_ids = pd.Index(points_df[instance_key].astype(str))
        if point_instance_ids.has_duplicates:
          raise ValueError(f"Points instance IDs for region '{region}' are not unique")
        candidate_ids.append((f"Points column '{instance_key}'", point_instance_ids, table_instance_ids))
    table_obs_ids = pd.Index(table_obs.index.astype(str))
    candidate_ids.append(("Points index", point_index, table_obs_ids))
    if not any(set(point_ids) == set(table_ids) for _, point_ids, table_ids in candidate_ids):
      raise ValueError(
        f"Points element '{region}' cannot be aligned to table rows by instance key or obs ID"
      )

    x_col, y_col = _infer_xy_columns(points_df)
    coordinates = points_df[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").to_numpy()
    if not np.isfinite(coordinates).all():
      raise ValueError(f"Points element '{region}' contains missing or non-finite coordinates")
    logger.info(
      f"Validated Points region '{region}': {len(points_df)} unique IDs with complete {x_col}/{y_col} coordinates"
    )

  return attrs, regions


def _restore_table_identity(adata, original_obs_names, original_uns):
  """Restore the source table order and metadata after coordinate conversion/concat."""
  adata.obs_names = adata.obs_names.astype(str)
  original_obs_names = pd.Index(original_obs_names.astype(str))
  if not original_obs_names.is_unique:
    raise ValueError("Original SpatialData table obs IDs are not unique")
  if not adata.obs_names.is_unique:
    raise ValueError(
      "Converted observations are not unique; coordinate-system subsets overlap or IDs were duplicated"
    )
  missing = original_obs_names.difference(adata.obs_names)
  extra = adata.obs_names.difference(original_obs_names)
  if len(adata) != len(original_obs_names) or len(missing) or len(extra):
    raise ValueError(
      "Converted table rows do not match the original SpatialData table "
      f"(original={len(original_obs_names)}, converted={len(adata)}, missing={len(missing)}, extra={len(extra)})"
    )
  adata = adata[original_obs_names].copy()
  for key, value in original_uns.items():
    adata.uns[key] = deepcopy(value)
  if not adata.obs_names.equals(original_obs_names):
    raise RuntimeError("Failed to restore the original SpatialData table observation order")
  if "spatial" not in adata.obsm or adata.obsm["spatial"].shape[0] != adata.n_obs:
    raise ValueError("Spatial conversion did not produce one coordinate row per table observation")
  spatial = np.asarray(adata.obsm["spatial"])
  if spatial.ndim != 2 or spatial.shape[1] < 2 or not np.isfinite(spatial[:, :2]).all():
    raise ValueError("adata.obsm['spatial'] contains incomplete or non-finite x/y coordinates")
  return adata


def _safe_condition_name(value):
  return str(value).replace(os.sep, "_").replace("/", "_").replace("\\", "_")

if args.cellcharter_col is None:
  args.cellcharter_col = "spatial_cluster"
if args.max_cluster is None:
  args.max_cluster = 10
if args.max_cluster < 2:
  raise ValueError("--max_cluster must be at least 2 for CellCharter AutoK")
if args.significance is None:
  args.significance = 0.05

concatenated_sdata = None
table_key = None
images = []
shapes = []
systems = []
if args.input_dir.endswith(".h5ad"):
  raise ValueError("CellCharter advance_analysis requires a SpatialData .zarr object; .h5ad input is not supported.")

log_step(logger, 1, 7, "loading SpatialData input")
concatenated_sdata = spd.read_zarr(args.input_dir)
table_key = next(iter(concatenated_sdata.tables.keys()))
original_table = concatenated_sdata.tables[table_key]
original_obs_names = pd.Index(original_table.obs_names.astype(str))
original_uns = deepcopy(original_table.uns)
spatialdata_attrs, annotated_regions = _validate_spatial_annotations(concatenated_sdata, original_table)
image_elements = list(concatenated_sdata.images.keys())
shape_elements = list(concatenated_sdata.shapes.keys())
valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
if len(valid_coord_systems) == 0:
  raise ValueError("SpatialData input has no coordinate system for spatial neighbor construction")
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
    data = to_legacy_anndata(
      concatenated_sdata,
      include_images=False,
      coordinate_system=s1,
      table_name=table_key,
    )
    if args.sample_col and args.sample_col not in data.obs.columns:
      data.obs[args.sample_col] = s1
    data.obsm['spatial'][:, 1] = np.max(data.obsm['spatial'][:, 1]) - data.obsm['spatial'][:, 1]
    adata_list.append(data)
  adata = ad.concat(adata_list, axis=0, merge='same', pairwise=True, index_unique=None)
else:
  coord_system = systems[0] if len(systems) > 0 else None
  point_annotated = any(region in concatenated_sdata.points for region in annotated_regions)
  adata = to_legacy_anndata(
    concatenated_sdata,
    include_images=not point_annotated,
    coordinate_system=coord_system,
    table_name=table_key,
  )
adata = _restore_table_identity(adata, original_obs_names, original_uns)
logger.info(f"Loaded {adata.n_obs} observations and {adata.n_vars} genes")

log_step(logger, 2, 7, "preparing counts layer")
counts_source, counts_source_name = _select_counts(adata)
x_contains_counts = _matrix_is_nonnegative_integer(adata.X)
adata.layers["counts"] = counts_source
if x_contains_counts:
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
logger.info(
  f"Using {counts_source_name} as scVI counts; "
  f"X normalization applied={x_contains_counts}"
)

if channel == "compare_analysis":
  if not args.sample_col:
    raise ValueError("--sample_col is required for compare_analysis")
  if args.sample_col not in adata.obs.columns and "library_id" in adata.obs.columns:
    adata.obs[args.sample_col] = adata.obs["library_id"]
  if args.sample_col not in adata.obs.columns:
    raise ValueError(
      f"Sample column '{args.sample_col}' is absent from the integrated table and no library_id fallback exists"
    )
  if adata.obs[args.sample_col].isna().any():
    raise ValueError(f"Sample column '{args.sample_col}' contains missing values")
  adata.obs[args.sample_col] = pd.Categorical(adata.obs[args.sample_col])
  adata.uns["spatial"] = {str(s): {} for s in adata.obs[args.sample_col].cat.categories}

if channel == "compare_analysis":
  scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.sample_col)
else:
  scvi.model.SCVI.setup_anndata(adata, layer="counts")


log_step(logger, 3, 7, "training SCVI model")
model = scvi.model.SCVI(adata)
model.train(early_stopping=True, enable_progress_bar=True)
adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

if channel == "compare_analysis":
  log_step(logger, 4, 7, "building spatial neighbors for compare analysis")
  sq.gr.spatial_neighbors(adata, library_key=args.sample_col, coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)
  spatial_graph = adata.obsp["spatial_connectivities"].tocoo()
  sample_codes = adata.obs[args.sample_col].cat.codes.to_numpy()
  if np.any(sample_codes[spatial_graph.row] != sample_codes[spatial_graph.col]):
    raise RuntimeError("Spatial neighbor graph contains edges between different samples")
  logger.info("Validated spatial neighbor graph: no cross-sample edges")
  cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter', sample_key=args.sample_col)
else:
  log_step(logger, 4, 7, "building spatial neighbors")
  sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)
  cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter')

cluster_trainer_params = {
  "accelerator": "gpu" if torch.cuda.is_available() else "cpu",
  "devices": 1,
  "strategy": "auto",
  "enable_progress_bar": False,
}

autok = cc.tl.ClusterAutoK(
  n_clusters=(2, args.max_cluster), 
  max_runs=10,
  convergence_tol=0.001,
  model_params={"trainer_params": cluster_trainer_params, "random_state": RANDOM_SEED},
)
log_step(logger, 5, 7, "running CellCharter spatial clustering")
autok.fit(adata, use_rep='X_cellcharter')
cc.pl.autok_stability(autok)
plt.savefig(
  os.path.join(dir_path, f"{sample_id}_celchar.png"),
  dpi=300,
  bbox_inches='tight')
plt.close()
best_k = int(autok.best_k)
predicted_labels = np.asarray(autok.predict(adata, use_rep='X_cellcharter')).astype(str)
adata.obs[args.cellcharter_col] = pd.Categorical(predicted_labels)
analysis_parameters = {
  "best_k": best_k,
  "random_seed": RANDOM_SEED,
  "counts_source": counts_source_name,
  "scvi_batch_key": args.sample_col if channel == "compare_analysis" else "",
  "spatial_coord_type": "generic",
  "spatial_delaunay": True,
  "spatial_percentile": 99,
  "spatial_library_key": args.sample_col if channel == "compare_analysis" else "",
  "aggregate_n_layers": 3,
  "aggregate_use_rep": "X_scVI",
  "aggregate_sample_key": args.sample_col if channel == "compare_analysis" else "",
  "autok_min_cluster": 2,
  "autok_max_cluster": args.max_cluster,
  "autok_max_runs": 10,
  "autok_convergence_tol": 0.001,
  "autok_accelerator": cluster_trainer_params["accelerator"],
  "threads": args.threads,
  "cellcharter_col": args.cellcharter_col,
}
adata.uns["cellcharter"] = analysis_parameters
logger.info(f"CellCharter selected best_k={best_k}; analysis parameters: {analysis_parameters}")
if concatenated_sdata is not None and table_key is not None:
  concatenated_sdata[table_key] = adata

if args.celltype_col and args.celltype_col in adata.obs.columns:
  logger.info(f"Running cell type enrichment with column {args.celltype_col}")
  cc.gr.enrichment(adata, group_key=args.cellcharter_col, label_key=args.celltype_col)
  cc.pl.enrichment(adata, group_key=args.cellcharter_col, label_key=args.celltype_col, figsize=(4,4), fontsize=6, dot_scale=1)
  plt.savefig(
    os.path.join(dir_path, f"{sample_id}_enrichment.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

condition_groups = []
replicate_counts = pd.Series(dtype=int)
if channel == "compare_analysis" and args.condition_col and args.condition_col not in adata.obs.columns:
  if "group" in adata.obs.columns:
    logger.warning(
      f"Condition column '{args.condition_col}' is absent; using project-standard column 'group' instead"
    )
    args.condition_col = "group"
  else:
    logger.warning(
      f"Condition column '{args.condition_col}' is absent and no 'group' fallback exists; "
      "skipping condition comparison"
    )

if channel == "compare_analysis" and args.condition_col and args.condition_col in adata.obs.columns:
  if adata.obs[args.condition_col].isna().any():
    raise ValueError(f"Condition column '{args.condition_col}' contains missing values")
  adata.obs[args.condition_col] = pd.Categorical(adata.obs[args.condition_col].astype(str))
  sample_condition = adata.obs[[args.sample_col, args.condition_col]].drop_duplicates()
  conditions_per_sample = sample_condition.groupby(args.sample_col, observed=True)[args.condition_col].nunique()
  if (conditions_per_sample > 1).any():
    invalid_samples = list(conditions_per_sample.index[conditions_per_sample > 1].astype(str))
    raise ValueError(
      "Each sample must map to exactly one condition; conflicting samples: "
      + ", ".join(invalid_samples)
    )
  condition_groups = list(adata.obs[args.condition_col].cat.categories)
  replicate_counts = sample_condition.groupby(args.condition_col, observed=True)[args.sample_col].nunique()
  logger.info(f"Condition-level biological replicate counts: {replicate_counts.to_dict()}")
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
      os.path.join(dir_path, f"{sample_id}_{_safe_condition_name(cond)}_enrichment.png"),
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
                        save_path=os.path.join(
                          dir_path,
                          f"{sample_id}_{_safe_condition_name(cond)}_Clusters_proportion.png"
                        ),
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

if channel == "compare_analysis" and len(condition_groups) >= 2:
  sufficient_replicates = bool((replicate_counts >= 2).all())
  adata.uns["cellcharter"]["condition_col"] = args.condition_col
  adata.uns["cellcharter"]["condition_groups"] = [str(value) for value in condition_groups]
  adata.uns["cellcharter"]["condition_replicates"] = {
    str(key): int(value) for key, value in replicate_counts.items()
  }
  adata.uns["cellcharter"]["diff_nhood_pvalues"] = sufficient_replicates
  adata.uns["cellcharter"]["diff_nhood_n_perms"] = 1000
  if sufficient_replicates:
    logger.info("Running sample-level differential neighborhood enrichment with 1,000 permutations")
  else:
    logger.warning(
      "At least one condition has fewer than two biological replicates; "
      "computing descriptive differential enrichment without significance testing"
    )
  cc.gr.diff_nhood_enrichment(
    adata,
    cluster_key=args.cellcharter_col,
    condition_key=args.condition_col,
    library_key=args.sample_col,
    pvalues=sufficient_replicates,
    n_jobs=args.threads,
    n_perms=1000)
  cc.pl.diff_nhood_enrichment(
    adata,
    cluster_key=args.cellcharter_col,
    condition_key=args.condition_col,
    condition_groups=condition_groups,
    annotate=True,
    figsize=(3,3),
    significance=args.significance if sufficient_replicates else None,
    fontsize=5)
  plt.savefig(
    os.path.join(dir_path, f"{sample_id}_diff_enrichment.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
elif channel == "compare_analysis" and len(condition_groups) == 1:
  logger.warning("Only one condition is present; skipping differential neighborhood enrichment")

def select_elements(items, keyword):
    if not keyword:
        return list(items)
    selected = [item for item in items if keyword in item]
    return selected if selected else list(items)
sample_cnt=adata.obs["sample"].astype(str).nunique() if "sample" in adata.obs else 1
concatenated_sdata[table_key] = adata
log_step(logger, 6, 7, "rendering CellCharter spatial plots")
active_vis_mode = "shape" if len(getattr(concatenated_sdata, "shapes", {})) > 0 else "point"
image_elements = list(concatenated_sdata.images.keys())
shape_elements = list(concatenated_sdata.shapes.keys())
point_elements = list(concatenated_sdata.points.keys())
valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
image_counts=len(image_elements)//sample_cnt
shape_count=len(shape_elements)//sample_cnt
shapes=shape_elements
images=[]
systems=[]
images = select_elements(image_elements, "hires")
systems = select_elements(valid_coord_systems, "hires")
region_name = adata.uns.get("spatialdata_attrs", {}).get("region")
if isinstance(region_name, (list, tuple)):
  region_name = region_name[0] if region_name else None
active_vis_mode, visual_region_name = resolve_visual_target(concatenated_sdata, region_name, "auto")
point_visual_elements = []
if active_vis_mode == "point":
    point_visual_elements = select_elements(point_elements, region_name)
    if not point_visual_elements:
        point_visual_elements = [visual_region_name] if visual_region_name in getattr(concatenated_sdata, "points", {}) else point_elements
    keep_ids = adata.obs_names.astype(str)
    point_transformations = {coord_system: Identity() for coord_system in valid_coord_systems} if valid_coord_systems else {"global": Identity()}
    for current_visual_region in point_visual_elements:
      source_points = concatenated_sdata.points[current_visual_region]
      points_df = source_points.compute() if hasattr(source_points, "compute") else source_points.copy()
      points_df = points_df.copy()
      points_df.index = points_df.index.astype(str)
      points_df = points_df.reindex(keep_ids).copy()
      points_df[args.cellcharter_col] = adata.obs[args.cellcharter_col].astype(str).reindex(points_df.index).fillna("unassigned").astype("category")
      points_df = points_df.dropna(subset=["x", "y"])
      concatenated_sdata.points[current_visual_region] = PointsModel.parse(
        points_df.sort_index(),
        coordinates={"x": "x", "y": "y"},
        transformations=point_transformations,
      )

for i in range(sample_cnt):
    sys = systems[0] if len(systems) < 2 else systems[i]
    title = images[i].replace("_hires_image", "")
    current_visual_region = visual_region_name
    if active_vis_mode == "point" and len(point_visual_elements) > 0:
      current_visual_region = point_visual_elements[0] if len(point_visual_elements) < 2 else point_visual_elements[min(i, len(point_visual_elements) - 1)]
    fig, axes = plt.subplots(2, 1, figsize=(20, 13))
    axes = np.ravel(axes)
    concatenated_sdata.pl.render_images(images[i]).pl.show(
      ax=axes[0],
      title="image",
      coordinate_systems=sys,
    )
    if active_vis_mode == "shape":
      concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(
        shapes[i],
        color=args.cellcharter_col,
        table_name=table_key,
        outline_alpha=0,
        fill_alpha=1,
      ).pl.show(
        ax=axes[1],
        coordinate_systems=sys,
        title=title,
      )
    else:
      concatenated_sdata.pl.render_images(images[i]).pl.render_points(
        current_visual_region,
        color=args.cellcharter_col,
        size=2.5,
        alpha=0.8,
      ).pl.show(
        ax=axes[1],
        coordinate_systems=sys,
        title=title,
      )
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

log_step(logger, 7, 7, f"saving CellCharter output to {args.output_zarr_path}")
concatenated_sdata.write(args.output_zarr_path)
logger.info("CellCharter module completed")
