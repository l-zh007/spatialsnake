###################preprocessing###########################################
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import numpy as np
import matplotlib.pyplot as plt
import argparse
import bbknn
import geosketch as sketch
from spatialsnake.workflow.function.pca_choosing import select_pca_dimensions
from spatialsnake.workflow.function.stereoseq_selection import resolve_stereoseq_table_key
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

logger = setup_logger("preprocess")

STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}


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
  if value_str in ["true", "1", "yes", "y", "t"]:
    return True
  if value_str in ["false", "0", "no", "n", "f", "none", "null", ""]:
    return False
  raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

parser = argparse.ArgumentParser(description='Filter, normalize, and prepare SpatialData for clustering')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Input SpatialData Zarr path')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Output preprocessed SpatialData Zarr path')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Sample or integrated-object identifier')

parser.add_argument('--type', type=str, required=True,
                   help='Spatial transcriptomics platform type')
parser.add_argument('--input_spec', type=str, required=False,
                   help='Stereo-seq input mode passed from sample.txt')

parser.add_argument('--min_counts', type=int, required=False,
                   help='Minimum UMI/count threshold per spot or cell')
parser.add_argument('--min_cells', type=int, required=False,
                   help='Minimum observations expressing a gene')
parser.add_argument('--variable', type=parse_bool, required=True,
                   help='Whether to select highly variable genes for PCA')
parser.add_argument('--filter_dict',type=str, required=False,
                   help='JSON mapping of sample-specific QC thresholds')
parser.add_argument('--NEIGHBORS', type=int, required=False,
                   help='Reserved neighbor setting from the preprocessing configuration')
parser.add_argument('--batch_method', type=str, required=False,
                   help='Batch correction method: None, harmony, or bbknn')
parser.add_argument('--mt_threshold', type=float, required=False,
                   help='Maximum mitochondrial-count percentage')
parser.add_argument('--n_top_genes', type=int, required=False,
                   help='Number of highly variable genes')
parser.add_argument('--n_comps', type=int, required=False,
                   help='Number of principal components to compute')
parser.add_argument('--sketch', type=parse_bool, required=False,
                   help='Whether to create a representative sketch for clustering')
parser.add_argument('--sample_rate', type=float, required=False,
                   help='Fraction retained within each stratum during sketching')
parser.add_argument('--threads', type=int, default=8,
                   help='Threads allocated to preprocessing')
args = parser.parse_args()
type=args.type
dir_path=os.path.dirname(args.output_zarr_path)
if args.threads < 1:
  parser.error("--threads must be >= 1")
sc.settings.n_jobs = args.threads

log_step(logger, 1, 7, "loading input data")
concatenated_sdata = spd.read_zarr(args.input_dir)
table = pick_table_key(concatenated_sdata, run_type=type, input_spec=args.input_spec)
adata = concatenated_sdata[table]
sample_count = adata.obs["sample"].astype(str).nunique() if "sample" in adata.obs else 1
logger.info(f"Loaded {adata.n_obs} observations and {adata.n_vars} genes from {sample_count} sample(s)")

log_step(logger, 2, 7, "calculating quality-control metrics")
# Calculate QC metrics for all cases
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"], 
    percent_top=(10, 20, 50),
    inplace=True, 
    log1p=True
)


if args.filter_dict is not None:
    log_step(logger, 3, 7, "filtering each sample with sample-specific thresholds")
    filter_dict = json.loads(args.filter_dict)
    candidate_columns = [column for column in ["sample", "region"] if column in adata.obs.columns]
    sample_column = None
    for column in candidate_columns:
      observed = set(adata.obs[column].astype(str).unique())
      if observed and observed.issubset(set(filter_dict)):
        sample_column = column
        break
    if sample_column is None:
      available = {
        column: sorted(adata.obs[column].astype(str).unique().tolist())
        for column in candidate_columns
      }
      raise RuntimeError(
        "Unable to match sample.txt filter IDs to an AnnData sample column. "
        f"filter IDs={sorted(filter_dict)}, observed={available}"
      )

    sample_labels = adata.obs[sample_column].astype(str)
    ordered_samples = list(dict.fromkeys(sample_labels.tolist()))
    unused_filters = sorted(set(filter_dict) - set(ordered_samples))
    if unused_filters:
      logger.warning(f"Ignoring filter thresholds for samples absent from the merged table: {','.join(unused_filters)}")

    position_col = "__spatialsnake_qc_position__"
    while position_col in adata.obs.columns:
      position_col = f"_{position_col}"
    adata.obs[position_col] = np.arange(adata.n_obs, dtype=int)
    kept_positions = []
    retained_gene_sets = []

    for sample_name in ordered_samples:
      if sample_name not in filter_dict:
        raise RuntimeError(f"Missing filter params for sample '{sample_name}' in filter_dict")
      sample_adata = adata[sample_labels == sample_name, :].copy()
      before_obs, before_vars = sample_adata.shape
      min_cells = int(filter_dict[sample_name][0])
      min_counts = int(filter_dict[sample_name][1])
      mt_threshold = float(filter_dict[sample_name][2])

      sc.pp.filter_cells(sample_adata, min_counts=min_counts)
      sample_adata = sample_adata[sample_adata.obs.pct_counts_mt < mt_threshold, :].copy()
      if sample_adata.n_obs == 0:
        raise RuntimeError(
          f"No observations remain for sample '{sample_name}' after min_counts={min_counts} "
          f"and mt_threshold={mt_threshold} filtering"
        )
      sc.pp.filter_genes(sample_adata, min_cells=min_cells)
      if sample_adata.n_vars == 0:
        raise RuntimeError(
          f"No genes remain for sample '{sample_name}' after min_cells={min_cells} filtering"
        )

      kept_positions.extend(sample_adata.obs[position_col].astype(int).tolist())
      retained_gene_sets.append(set(sample_adata.var_names.astype(str)))
      logger.info(
        f"Sample {sample_name}: {before_obs}->{sample_adata.n_obs} observations, "
        f"{before_vars}->{sample_adata.n_vars} genes "
        f"(min_counts={min_counts}, min_cells={min_cells}, mt_threshold={mt_threshold})"
      )

    if not kept_positions:
      raise RuntimeError("No observations remain after per-sample filtering")
    common_genes = set.intersection(*retained_gene_sets)
    if not common_genes:
      raise RuntimeError(
        "No common genes remain after applying the per-sample min_cells thresholds; "
        "reduce min_cells or inspect sample quality"
      )

    # Subsetting the original merged AnnData is equivalent to an inner-join
    # concatenation of the independently filtered samples, while preserving
    # SpatialData table metadata, layers, obsm and geometry instance links.
    keep_position_set = set(kept_positions)
    keep_obs = np.fromiter(
      (position in keep_position_set for position in range(adata.n_obs)),
      dtype=bool,
      count=adata.n_obs,
    )
    keep_vars = np.fromiter(
      (str(gene) in common_genes for gene in adata.var_names),
      dtype=bool,
      count=adata.n_vars,
    )
    adata = adata[keep_obs, keep_vars].copy()
    del adata.obs[position_col]
    logger.info(
      f"Merged {len(ordered_samples)} independently filtered samples using "
      f"{adata.n_vars} genes retained in every sample"
    )
else:
    min_counts=args.min_counts
    min_cells=args.min_cells
    log_step(logger, 3, 7, f"filtering cells and genes with min_counts={min_counts}, min_cells={min_cells}")
    sc.pp.filter_cells(adata,min_counts=min_counts)
    maxi=float(args.mt_threshold)
    adata=adata[adata.obs.pct_counts_mt<maxi, :].copy()
    sc.pp.filter_genes(adata,min_cells=min_cells)

logger.info(f"After filtering: {adata.n_obs} observations and {adata.n_vars} genes")

log_step(logger, 4, 7, "saving quality-control plots")
sc.pl.violin(adata=adata, keys=["log1p_total_counts"], stripplot=False, inner="box",show=False, groupby="region",)
plt.title("Total UMI by Sample")
if args.filter_dict is None and args.min_counts is not None:
  plt.axhline(
    y=np.log1p(args.min_counts),
    color="r",
    linestyle="--",
    linewidth=1,
    label=f"min_counts={args.min_counts}",
  )
  plt.legend(frameon=False, fontsize=8)
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}filtered_Total_UMI.png"),
  dpi=300,
  bbox_inches='tight')
plt.close()

sc.pl.violin(adata=adata, keys=["log1p_n_genes_by_counts"], groupby="region", stripplot=False, inner="box",show=False)
plt.title("Total Genes by Sample")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}filtered_Total_Genes.png"),
    dpi=300,
    bbox_inches='tight')
plt.close()

sc.pl.violin(adata=adata, keys=["log1p_total_counts_mt"], groupby="region", stripplot=False, inner="box",show=False)
plt.title("Mitochondrial Genes by Sample")
plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}_Mitochondrial_Genes.png"),
  dpi=300,
  bbox_inches='tight')
plt.close()


sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_scatter.png"),
    dpi=300,
    bbox_inches='tight')
plt.close()


if "counts" not in adata.layers:
    adata.raw = adata.copy()
else:
    raw_adata = adata.copy()
    raw_adata.X = adata.layers['counts']
    adata.raw = raw_adata
    del raw_adata

log_step(logger, 5, 7, "normalizing expression and running PCA")
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)



use_hvg = False
if args.variable:
  sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,n_top_genes=args.n_top_genes)
  sc.pl.highly_variable_genes(adata, show=False)
  plt.savefig(
  os.path.join(dir_path, f"{args.sample_id}_highly_variable.png"),
  dpi=300,
  bbox_inches='tight')
  plt.close()
  use_hvg = True
  logger.info(f"Selected highly variable genes with n_top_genes={args.n_top_genes}")

# sc.pp.scale is intentionally omitted because the existing workflow performs
# PCA on the log-normalized expression matrix.
pca_feature_count = int(adata.var["highly_variable"].sum()) if use_hvg else adata.n_vars
max_pca_components = min(adata.n_obs, pca_feature_count) - 1
if max_pca_components < 1:
  raise RuntimeError(
    "PCA requires at least two observations and two selected genes after filtering"
  )
requested_n_comps = int(50 if args.n_comps is None else args.n_comps)
if requested_n_comps < 1:
  raise ValueError("--n_comps must be >= 1")
effective_n_comps = min(requested_n_comps, max_pca_components)
if effective_n_comps != requested_n_comps:
  logger.warning(
    f"Requested n_comps={requested_n_comps}, but the filtered data support at most "
    f"{max_pca_components}; using n_comps={effective_n_comps}"
  )
sc.tl.pca(
  adata,
  n_comps=effective_n_comps,
  use_highly_variable=use_hvg,
  random_state=0,
)

log_step(logger, 6, 7, f"handling batch effect with method={args.batch_method}")
batch_key = next(
  (key for key in ("sample", "region") if key in adata.obs and adata.obs[key].nunique() > 1),
  None,
)
requested_batch_method = str(args.batch_method or "None").strip().lower()
if requested_batch_method in {"", "none", "false", "null"}:
  requested_batch_method = "none"
if requested_batch_method not in {"none", "harmony", "bbknn"}:
  raise ValueError(
    f"Unsupported batch_method '{args.batch_method}'. Use None, harmony, or bbknn."
  )

applied_batch_method = requested_batch_method
if batch_key is None:
  applied_batch_method = "none"
  logger.info("Batch correction skipped because fewer than two batches were detected")
elif requested_batch_method == "none":
  logger.info("Batch correction disabled")
elif requested_batch_method == "harmony":
  logger.info(f"Running Harmony integration with batch_key={batch_key}")
  sce.pp.harmony_integrate(
    adata,
    key=batch_key,
    basis="X_pca",
    adjusted_basis="X_pca_harmony",
    max_iter_harmony=20,
  )

run_sketch = args.sketch
if run_sketch:
  if args.sample_rate is None or not 0 < args.sample_rate <= 1:
    raise ValueError("--sample_rate must be in the interval (0, 1]")
  sketch_key = batch_key
  if sketch_key is None:
    sketch_key = next((key for key in ("sample", "region", "group") if key in adata.obs), None)
  if sketch_key is None:
    sketch_labels = np.repeat("all", adata.n_obs)
  else:
    sketch_labels = adata.obs[sketch_key].astype(str).to_numpy()
  sketch_rep = "X_pca_harmony" if applied_batch_method == "harmony" else "X_pca"
  logger.info(
    f"Running sketch sampling with sample_rate={args.sample_rate}, "
    f"stratification={sketch_key or 'none'}, representation={sketch_rep}"
  )
  selected_positions = []
  for value in dict.fromkeys(sketch_labels.tolist()):
    positions = np.flatnonzero(sketch_labels == value)
    target_n = max(1, min(len(positions), int(round(args.sample_rate * len(positions)))))
    if target_n == len(positions):
      local_index = np.arange(len(positions), dtype=int)
    else:
      local_index = sketch.gs(
        X=adata.obsm[sketch_rep][positions],
        N=target_n,
        seed=1,
      )
    selected_positions.extend(positions[np.asarray(local_index, dtype=int)].tolist())
  selected_positions = np.asarray(sorted(set(selected_positions)), dtype=int)
  sdata = adata[selected_positions, :].copy()
  logger.info(f"Sketch retained {sdata.n_obs} of {adata.n_obs} observations")
else:
  sdata = adata

if applied_batch_method == "bbknn":
  logger.info(f"Running BBKNN integration with batch_key={batch_key}")
  bbknn.bbknn(
    sdata,
    batch_key=batch_key,
    use_rep="X_pca",
    n_pcs=effective_n_comps,
  )

preprocess_metadata = {
  "batch_method_requested": requested_batch_method,
  "batch_method": applied_batch_method,
  "batch_key": batch_key or "",
  "clustering_rep": "X_pca_harmony" if applied_batch_method == "harmony" else "X_pca",
  "use_precomputed_neighbors": applied_batch_method == "bbknn",
  "sketch": bool(run_sketch),
  "full_n_obs": int(adata.n_obs),
  "sketch_n_obs": int(sdata.n_obs),
  "requested_n_comps": requested_n_comps,
  "effective_n_comps": effective_n_comps,
}
sdata.uns["spatialsnake_preprocess"] = preprocess_metadata
if run_sketch:
  adata.uns["spatialsnake_preprocess"] = preprocess_metadata
  adata.write(os.path.join(dir_path, "sketch.h5ad"))
  logger.info("Saved the full preprocessed query to sketch.h5ad")
sc.pl.pca_variance_ratio(sdata, log=True, n_pcs=min(20, effective_n_comps), show=False)
plt.title("pca_variance_ratio")
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}pca_variance_ratio.png"),
    dpi=300,
    bbox_inches='tight')
plt.close()

pca_selection=select_pca_dimensions(sdata)
logger.info(f"Recommended PCs: {pca_selection}")

log_step(logger, 7, 7, f"saving preprocessed data to {args.output_zarr_path}")
spatial_attrs = sdata.uns.get('spatialdata_attrs', {})
instance_key = spatial_attrs.get('instance_key')
sdata.obs[instance_key] = sdata.obs[instance_key].astype(str)
sdata.obs['region'] = sdata.obs['region'].astype('category')
concatenated_sdata[table]=sdata
concatenated_sdata.write(os.path.join(args.output_zarr_path),overwrite=True)
logger.info("Preprocess module completed")
