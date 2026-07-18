import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as spd
from sklearn.cluster import KMeans

from spatialsnake.workflow.function.logging_utils import log_step, setup_logger
from spatialsnake.workflow.function.stereoseq_selection import (
  extract_table_spatial_metadata,
  resolve_stereoseq_table_key,
)


STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}
CLUSTER_REP_KEY = "X_spatialsnake_cluster"
logger = setup_logger("clustering")


def pick_table_key(sdata, run_type=None, input_spec=None):
  table_keys = list(getattr(sdata, "tables", {}).keys())
  if not table_keys:
    raise RuntimeError("No table found in SpatialData input.")
  if run_type in STEREOSEQ_TYPES:
    return resolve_stereoseq_table_key(input_spec, table_keys)
  return "table" if "table" in table_keys else table_keys[0]


def parse_bool(value):
  if isinstance(value, bool):
    return value
  if value is None:
    return False
  value_str = str(value).strip().lower()
  if value_str in {"true", "1", "yes", "y", "t"}:
    return True
  if value_str in {"false", "0", "no", "n", "f", "none", "null", ""}:
    return False
  raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def get_preprocess_metadata(adata):
  metadata = adata.uns.get("spatialsnake_preprocess", {})
  return dict(metadata) if hasattr(metadata, "items") else {}


def select_representation(adata, metadata):
  batch_method = str(metadata.get("batch_method", "none")).lower()
  representation = str(metadata.get("clustering_rep", "X_pca"))
  if not metadata and "X_pca_harmony" in adata.obsm:
    batch_method = "harmony"
    representation = "X_pca_harmony"
  if representation not in adata.obsm:
    raise RuntimeError(
      f"Clustering representation '{representation}' is missing from obsm. "
      "Rerun preprocessing so that batch correction and clustering use compatible data."
    )
  return batch_method, representation


def add_clustering_representation(adata, representation, requested_pcs):
  available_pcs = int(adata.obsm[representation].shape[1])
  if requested_pcs < 1:
    raise ValueError("--pcs must be >= 1")
  effective_pcs = min(requested_pcs, available_pcs)
  if effective_pcs != requested_pcs:
    logger.warning(
      f"Requested pcs={requested_pcs}, but representation '{representation}' contains "
      f"{available_pcs} components; using pcs={effective_pcs}"
    )
  adata.obsm[CLUSTER_REP_KEY] = np.asarray(
    adata.obsm[representation][:, :effective_pcs], dtype=np.float32
  ).copy()
  return effective_pcs


parser = argparse.ArgumentParser(description="Cluster a preprocessed SpatialData object")
parser.add_argument("--input_dir", type=str, required=True, help="Input preprocessed SpatialData Zarr path")
parser.add_argument("--output_zarr_path", type=str, required=True, help="Output clustered SpatialData Zarr path")
parser.add_argument("--sample_id", type=str, required=True, help="Sample or integrated-object identifier")
parser.add_argument("--type", type=str, required=True, help="Spatial transcriptomics platform type")
parser.add_argument("--tsne", type=parse_bool, required=False, help="Whether to compute a t-SNE embedding")
parser.add_argument("--MIN_DIST", type=float, required=False, help="UMAP minimum distance")
parser.add_argument("--SPREAD", type=float, required=False, help="UMAP spread")
parser.add_argument("--RES", type=float, required=False, help="Leiden or Louvain resolution")
parser.add_argument(
  "--cluster_algorithm",
  type=str,
  required=False,
  help="Clustering algorithm: leiden, louvain, or kmeans",
)
parser.add_argument("--n_clusters", type=float, required=False, help="Number of clusters for KMeans")
parser.add_argument("--sketch", type=parse_bool, required=False, help="Whether preprocessing used sketch sampling")
parser.add_argument("--pcs", type=int, required=False, help="Number of components used for clustering")
parser.add_argument("--NEIGHBORS", type=int, required=False, help="Number of neighbors in the graph")
parser.add_argument("--input_spec", type=str, required=False, help="Stereo-seq input mode passed from sample.txt")
parser.add_argument("--threads", type=int, default=8, help="Threads allocated to clustering")
args = parser.parse_args()

run_type = args.type
cluster_algorithm = (args.cluster_algorithm or "leiden").lower()
if args.threads < 1:
  parser.error("--threads must be >= 1")
if args.NEIGHBORS is None or args.NEIGHBORS < 2:
  parser.error("--NEIGHBORS must be >= 2")
sc.settings.n_jobs = args.threads

dir_path = os.path.dirname(args.output_zarr_path)

log_step(logger, 1, 6, "loading preprocessed data")
concatenated_sdata = spd.read_zarr(args.input_dir)
table = pick_table_key(concatenated_sdata, run_type=run_type, input_spec=args.input_spec)
adata_for_sketch = concatenated_sdata[table]
metadata = get_preprocess_metadata(adata_for_sketch)
batch_method, clustering_rep = select_representation(adata_for_sketch, metadata)
requested_pcs = int(20 if args.pcs is None else args.pcs)
effective_pcs = add_clustering_representation(
  adata_for_sketch, clustering_rep, requested_pcs
)

logger.info(
  f"Loaded {adata_for_sketch.n_obs} observations and {adata_for_sketch.n_vars} genes; "
  f"batch_method={batch_method}, representation={clustering_rep}, pcs={effective_pcs}"
)

log_step(logger, 2, 6, f"preparing the neighbor graph with {effective_pcs} components")
use_precomputed_neighbors = bool(metadata.get("use_precomputed_neighbors", False))
if batch_method == "bbknn" and use_precomputed_neighbors:
  if "neighbors" not in adata_for_sketch.uns or not {
    "connectivities",
    "distances",
  }.issubset(adata_for_sketch.obsp):
    raise RuntimeError(
      "Preprocessing metadata specifies BBKNN, but its neighbor graph is missing. "
      "Rerun preprocessing before clustering."
    )
  logger.info("Using the BBKNN graph generated during preprocessing")
else:
  sc.pp.neighbors(
    adata_for_sketch,
    n_neighbors=args.NEIGHBORS,
    use_rep=CLUSTER_REP_KEY,
    metric="euclidean",
    random_state=0,
  )

metadata.update(
  {
    "clustering_rep": clustering_rep,
    "requested_pcs": requested_pcs,
    "effective_pcs": effective_pcs,
    "neighbor_graph": "bbknn" if batch_method == "bbknn" and use_precomputed_neighbors else "scanpy",
  }
)
adata_for_sketch.uns["spatialsnake_preprocess"] = metadata

run_sketch = bool(metadata.get("sketch", args.sketch))
if args.tsne and run_sketch:
  logger.warning(
    "t-SNE is skipped in sketch mode because Scanpy ingest does not support "
    "projecting a reference t-SNE embedding to the full dataset"
  )
elif args.tsne:
  logger.info("Running t-SNE embedding")
  perplexity = min(30, max(1, (adata_for_sketch.n_obs - 1) // 3))
  sc.tl.tsne(
    adata_for_sketch,
    use_rep=CLUSTER_REP_KEY,
    perplexity=perplexity,
    random_state=0,
  )

log_step(logger, 3, 6, f"running {cluster_algorithm} clustering")
if cluster_algorithm == "leiden":
  sc.tl.leiden(
    adata_for_sketch,
    flavor="igraph",
    key_added="clusters",
    resolution=args.RES,
    random_state=0,
  )
elif cluster_algorithm == "louvain":
  sc.tl.louvain(
    adata_for_sketch,
    resolution=args.RES,
    key_added="clusters",
    random_state=0,
  )
elif cluster_algorithm == "kmeans":
  if args.n_clusters is None or int(args.n_clusters) < 2:
    raise ValueError("--n_clusters must be >= 2 for KMeans")
  kmeans = KMeans(n_clusters=int(args.n_clusters), random_state=0, n_init=10)
  labels = kmeans.fit_predict(adata_for_sketch.obsm[CLUSTER_REP_KEY])
  adata_for_sketch.obs["clusters"] = pd.Categorical(labels.astype(str))
else:
  raise ValueError(f"Unsupported cluster_algorithm: {args.cluster_algorithm}")

adata_for_sketch.obs["clusters"] = adata_for_sketch.obs["clusters"].astype("category")

log_step(logger, 4, 6, "running UMAP embedding")
sc.tl.umap(
  adata_for_sketch,
  min_dist=args.MIN_DIST,
  spread=args.SPREAD,
  random_state=0,
)

if run_sketch:
  logger.info("Projecting sketch clusters and UMAP coordinates back to the full data")
  full_path = os.path.join(os.path.dirname(args.input_dir), "sketch.h5ad")
  if not os.path.exists(full_path):
    raise FileNotFoundError(
      f"Sketch mode requires the full preprocessed query at {full_path}"
    )
  full_adata = sc.read_h5ad(full_path)
  if not full_adata.var_names.equals(adata_for_sketch.var_names):
    raise RuntimeError("The sketch reference and full query have different gene order")
  if clustering_rep not in full_adata.obsm:
    raise RuntimeError(
      f"The full sketch query does not contain '{clustering_rep}'. "
      "Rerun preprocessing with the updated batch-correction workflow."
    )
  if metadata.get("full_n_obs") and full_adata.n_obs != int(metadata["full_n_obs"]):
    raise RuntimeError(
      "The full sketch query observation count does not match preprocessing metadata"
    )

  full_adata.obsm[CLUSTER_REP_KEY] = np.asarray(
    full_adata.obsm[clustering_rep][:, :effective_pcs], dtype=np.float32
  ).copy()
  cluster_values = np.empty(full_adata.n_obs, dtype=object)
  umap_values = np.empty(
    (full_adata.n_obs, adata_for_sketch.obsm["X_umap"].shape[1]), dtype=np.float32
  )
  chunk_size = 100000
  for start in range(0, full_adata.n_obs, chunk_size):
    end = min(start + chunk_size, full_adata.n_obs)
    query_chunk = full_adata[start:end].copy()
    sc.tl.ingest(
      query_chunk,
      adata_for_sketch,
      obs="clusters",
      embedding_method=("umap",),
    )
    cluster_values[start:end] = query_chunk.obs["clusters"].astype(str).to_numpy()
    umap_values[start:end] = query_chunk.obsm["X_umap"]

  categories = adata_for_sketch.obs["clusters"].cat.categories.astype(str)
  full_adata.obs["clusters"] = pd.Categorical(cluster_values, categories=categories)
  full_adata.obsm["X_umap"] = umap_values
  if "umap" in adata_for_sketch.uns:
    full_adata.uns["umap"] = adata_for_sketch.uns["umap"]
  full_adata.uns["spatialsnake_preprocess"] = metadata
  del full_adata.obsm[CLUSTER_REP_KEY]
  adata_clustered = full_adata
else:
  adata_clustered = adata_for_sketch

if adata_clustered.obs["clusters"].isna().any():
  raise RuntimeError("Cluster transfer produced missing labels in the final table")
adata_clustered.obs["clusters"] = adata_clustered.obs["clusters"].astype("category")

log_step(logger, 5, 6, "saving clustering diagnostic plots")
sc.pl.umap(
  adata_clustered,
  color=["total_counts", "n_genes_by_counts", "clusters"],
  wspace=0.4,
  show=False,
)
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}UMAP.png"),
  dpi=300,
  bbox_inches="tight",
)
plt.close()

if batch_method == "harmony":
  batch_key = str(metadata.get("batch_key", ""))
  plot_key = None
  if "sample" in adata_clustered.obs and adata_clustered.obs["sample"].nunique() > 1:
    plot_key = "sample"
  elif batch_key in adata_clustered.obs and adata_clustered.obs[batch_key].nunique() > 1:
    plot_key = batch_key
  if plot_key is None:
    logger.warning("Harmony sample UMAP was skipped because no multi-level batch column is available")
  else:
    sc.pl.umap(
      adata_clustered,
      color=plot_key,
      title=f"Harmony-corrected UMAP by {plot_key}",
      frameon=False,
      legend_loc="right margin",
      show=False,
    )
    plt.savefig(
      os.path.join(dir_path, f"{args.sample_id}_harmony_sample_UMAP.png"),
      dpi=300,
      bbox_inches="tight",
      facecolor="white",
    )
    plt.close()
    logger.info(f"Saved the Harmony sample UMAP colored by {plot_key}")

if args.tsne and not run_sketch:
  sc.pl.tsne(
    adata_clustered,
    color=["total_counts", "n_genes_by_counts", "clusters"],
    ncols=3,
    show=False,
  )
  plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}tsne.png"),
    dpi=300,
    bbox_inches="tight",
  )
  plt.close()

distribution_key = next(
  (key for key in ("region", "sample") if key in adata_clustered.obs), None
)
if distribution_key is None:
  logger.warning("Cluster distribution plot skipped because region/sample metadata are absent")
else:
  distribution = pd.crosstab(
    adata_clustered.obs[distribution_key], adata_clustered.obs["clusters"]
  )
  figure_height = max(3, 0.35 * distribution.shape[0] + 1.5)
  figure_width = max(5, 0.35 * distribution.shape[1] + 2)
  fig, ax = plt.subplots(figsize=(figure_width, figure_height))
  image = ax.imshow(distribution.to_numpy(), cmap="hot", interpolation="nearest", aspect="auto")
  ax.set_title("Cell Distribution Across Clusters")
  ax.set_xlabel("Cluster")
  ax.set_ylabel(distribution_key)
  ax.set_xticks(range(distribution.shape[1]), distribution.columns.astype(str), rotation=90)
  ax.set_yticks(range(distribution.shape[0]), distribution.index.astype(str))
  fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
  fig.tight_layout()
  fig.savefig(
    os.path.join(dir_path, f"{args.sample_id}Cell_Distribution_Across_Clusters.png"),
    dpi=300,
    bbox_inches="tight",
  )
  plt.close(fig)

log_step(logger, 6, 6, f"saving clustered SpatialData to {args.output_zarr_path}")
table_meta = extract_table_spatial_metadata(adata_clustered)
instance_key = table_meta["instance_key"]
adata_clustered.obs_names = adata_clustered.obs_names.astype(str)
adata_clustered.obs[instance_key] = adata_clustered.obs[instance_key].astype(str)

if not run_sketch and batch_method == "bbknn" and use_precomputed_neighbors:
  # BBKNN uses X_pca directly; the temporary sliced representation is not
  # referenced by its stored graph.
  del adata_clustered.obsm[CLUSTER_REP_KEY]

concatenated_sdata[table] = adata_clustered
concatenated_sdata.write(args.output_zarr_path, overwrite=True)
logger.info(
  f"Clustering module completed with {adata_clustered.n_obs} observations and "
  f"{adata_clustered.obs['clusters'].nunique()} clusters"
)
