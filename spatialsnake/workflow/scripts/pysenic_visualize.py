"""Preserve the pySCENIC rule/API while isolating scoring and plotting memory.

The existing full and AUCell H5AD outputs are exact, reusable checkpoints.
Only an internal cache manifest is added; no scientific output is removed.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")
import argparse
import gc
import hashlib
import importlib.metadata
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True,
                   help='Path to the raw data directory')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--celltype', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--outputs', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--num_workers', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--types', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--regulons', type=str, required=True,
                   help='Path to the raw data directory')
parser.add_argument('--top_regulons', type=int, default=20,
                   help='Number of regulons shown in dotplot and violin plots. Use <=0 to show all regulons.')
parser.add_argument('--min_regulon_genes', type=int, default=10,
                   help='Minimum number of target genes required for a regulon to be used in AUCell.')
parser.add_argument('--stage', choices=['all', 'compute', 'plot', 'export'], default='all', help=argparse.SUPPRESS)
args = parser.parse_args()
output_dir = os.path.dirname(args.outputs)
os.makedirs(output_dir, exist_ok=True)
top_regulons = max(0, int(args.top_regulons))
min_regulon_genes = max(1, int(args.min_regulon_genes))
full_checkpoint = Path(output_dir) / f"{args.sample_id}.h5ad"
auc_checkpoint = Path(output_dir) / f"{args.sample_id}_aucell.h5ad"
cache_manifest = Path(output_dir) / f".{args.sample_id}.aucell_cache.json"


def figure_pages(filename, labels, per_page):
    path = Path(output_dir) / filename
    for old_page in path.parent.glob(f"{path.stem}_page[0-9]*{path.suffix}"):
        old_page.unlink()
    for start in range(0, len(labels), per_page):
        page = start // per_page + 1
        target = path if page == 1 else path.with_name(f"{path.stem}_page{page:02d}{path.suffix}")
        yield target, list(labels[start:start + per_page])


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_digest(path):
    path = Path(path)
    if path.is_file():
        return file_digest(path)
    if not path.is_dir():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    for entry in sorted(p for p in path.rglob("*") if p.is_file()):
        digest.update(str(entry.relative_to(path)).encode())
        digest.update(b"\0")
        digest.update(file_digest(entry).encode())
    return digest.hexdigest()


def compute_identity():
    versions = {}
    for package in ("pyscenic", "ctxcore", "numpy", "pandas", "scipy", "anndata", "scanpy", "spatialdata"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return {
        "schema": 1,
        "input_path": str(Path(args.input_dir).resolve()),
        "input_sha256": input_digest(args.input_dir),
        "regulons_sha256": file_digest(args.regulons),
        "min_regulon_genes": min_regulon_genes,
        "num_workers": args.num_workers,
        "source_sha256": file_digest(__file__),
        "compat_sha256": file_digest(Path(__file__).with_name("pyscenic_numpy_compat.py")),
        "versions": versions,
    }


def read_valid_cache():
    try:
        record = json.loads(cache_manifest.read_text())
        if record["identity"] != compute_identity():
            return None
        if any(file_digest(path) != record["checkpoints"][path.name]
               for path in (full_checkpoint, auc_checkpoint)):
            return None
        return record
    except (OSError, ValueError, KeyError, TypeError):
        return None


def run_stage(stage):
    subprocess.run([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], '--stage', stage], check=True)


if args.stage == 'all':
    if read_valid_cache() is None:
        run_stage('compute')
    else:
        print('Reusing validated AUCell checkpoints; inference is unchanged.', flush=True)
    # Same computation -> visualization -> final Loom order inside the same rule.
    run_stage('plot')
    run_stage('export')
    sys.exit(0)

# Heavy scientific imports are confined to the stage process, never the controller.
import spatialdata as spd
import scanpy as sc
import spatialdata_plot as splt
import pandas as pd
import numpy as np
import re, pickle, math
from functools import partial
from collections import OrderedDict
from cytoolz import compose
import anndata as ad
import h5py
from anndata.io import read_elem
from pyscenic_numpy_compat import apply_numpy_compat, apply_pkg_resources_compat

apply_numpy_compat()
apply_pkg_resources_compat()

from pyscenic.export import add_scenic_metadata
from pyscenic.utils import load_motifs
try:
  from pyscenic.utils import load_regulons
except Exception:
  load_regulons = None
from pyscenic.transform import df2regulons
from pyscenic.aucell import create_rankings, aucell4r, DTYPE
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from spatialsnake.workflow.function.plot import plot_auc_heatmap_scanpy
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step
import argparse
logger = setup_logger("pyscenic_visualize")

def read_regulons():
    if load_regulons is not None:
      regulons = load_regulons(args.regulons)
    else:
      motifs = load_motifs(args.regulons)
      if isinstance(motifs.columns, pd.MultiIndex):
        motifs.columns = motifs.columns.droplevel(0)
      regulons = list(filter(lambda r: len(r) >= min_regulon_genes, df2regulons(motifs[(motifs['NES'] >= 3.0)])))
      regulons = list(map(lambda r: r.rename(r.transcription_factor), regulons))
    regulons = [reg for reg in regulons if len(reg) >= min_regulon_genes]
    logger.info(f"Loaded {len(regulons)} regulon(s) with at least {min_regulon_genes} target genes")
    if len(regulons) == 0:
      raise ValueError(
          f"No regulons passed the minimum target-gene threshold ({min_regulon_genes}). "
          "Lower pyscenic_min_regulon_genes or check the ctx regulon output."
      )

    return regulons


def read_expression_table():
    if os.path.splitext(args.input_dir)[1].lower() == ".h5ad":
        return sc.read_h5ad(args.input_dir)
    sdata = spd.read_zarr(args.input_dir)
    return sdata.tables[list(sdata.tables)[-1]]


def rank_expression(expression, cells, genes):
    # Ranking is independent for each row. Reuse the same random gene shuffle
    # for every block, including ties, and advance the RNG only once as before.
    if expression.shape[0] == 0:
        values = expression.toarray() if issparse(expression) else expression
        return create_rankings(pd.DataFrame(values, index=cells, columns=genes))
    block_size = max(1, (64 * 1024**2) // max(1, expression.shape[1] * expression.dtype.itemsize))
    values = np.empty(expression.shape, dtype=DTYPE)
    state = np.random.get_state()
    columns = None
    for start in range(0, expression.shape[0], block_size):
        stop = min(start + block_size, expression.shape[0])
        block = expression[start:stop]
        if issparse(block):
            block = block.toarray()
        np.random.set_state(state)
        ranked = create_rankings(pd.DataFrame(block, index=cells[start:stop], columns=genes))
        if columns is None:
            columns = ranked.columns
            next_state = np.random.get_state()
        values[start:stop] = ranked.to_numpy()
        del block, ranked
    np.random.set_state(next_state)
    return pd.DataFrame(values, index=cells, columns=columns)


if args.stage == 'compute':
    cache_manifest.unlink(missing_ok=True)
    identity = compute_identity()
    log_step(logger, 1, 5, "loading expression data and regulons")
    adata = read_expression_table()

    regulons = read_regulons()
    log_step(logger, 2, 5, "calculating AUCell matrix")
    exp_matrix = adata.X
    gene_names = adata.var_names
    cell_names = adata.obs_names
    # Ranking only needs X and its axes, not raw/layers/images. Reload the same
    # table after scoring to preserve every field of the existing H5AD output.
    del adata
    gc.collect()
    if issparse(exp_matrix):
      nnz = exp_matrix.nnz
      dense_bytes = exp_matrix.shape[0] * exp_matrix.shape[1] * exp_matrix.dtype.itemsize
      logger.info(
          f"Sparse expression matrix detected with shape={exp_matrix.shape}, nnz={nnz}, "
          f"estimated dense {exp_matrix.dtype} memory={dense_bytes / 1024**3:.2f} GB"
      )
      # uint32 ranks retain the original float32-based size limit, including
      # for float64 inputs; only the reported expression estimate changes.
      if exp_matrix.shape[0] * exp_matrix.shape[1] * np.dtype(DTYPE).itemsize > 8 * 1024**3:
        raise MemoryError(
            "The AUCell ranking matrix is estimated to exceed the existing 8 GB size limit. "
            "Please reduce object size before visualization."
        )
    else:
      logger.info(f"Dense expression matrix detected with shape={exp_matrix.shape}")
    rankings = rank_expression(exp_matrix, cell_names, gene_names)
    del exp_matrix
    gc.collect()
    auc_mtx = aucell4r(rankings, regulons, num_workers=args.num_workers)
    del rankings
    gc.collect()

    adata = read_expression_table()
    add_scenic_metadata(adata, auc_mtx, regulons)
    adata.write(os.path.join(output_dir,f"{args.sample_id}.h5ad"))

    aucell_adata = sc.AnnData(X=auc_mtx.loc[adata.obs_names].values)
    aucell_adata.obs = adata.obs
    aucell_adata.var_names = auc_mtx.columns
    aucell_adata.write(os.path.join(output_dir,f"{args.sample_id}_aucell.h5ad"))


    rng = np.random.get_state()
    record = {
        "identity": identity,
        "checkpoints": {path.name: file_digest(path) for path in (full_checkpoint, auc_checkpoint)},
        "auc_index": [str(value) for value in auc_mtx.index],
        "auc_index_name": auc_mtx.index.name,
        "auc_columns_name": auc_mtx.columns.name,
        "numpy_random_state": [rng[0], rng[1].tolist(), rng[2], rng[3], rng[4]],
        "python_random_state": random.getstate(),
    }
    with tempfile.NamedTemporaryFile(mode='w', dir=output_dir, prefix='.aucell_cache-', suffix='.json', delete=False) as stream:
        json.dump(record, stream)
        temporary = Path(stream.name)
    os.replace(temporary, cache_manifest)
    logger.info("AUCell checkpoints saved; releasing scoring process memory")
    sys.exit(0)

if args.stage == 'export':
    # Loom requires expression, but no dense ranking matrix or plotting objects.
    adata = sc.read_h5ad(full_checkpoint)
    log_step(logger, 5, 5, f"writing pySCENIC loom output to {args.outputs}")
    adata.write_loom(args.outputs)
    logger.info("pySCENIC visualization module completed")
    sys.exit(0)

# Plot only activity scores and original observation metadata, not full expression.
record = json.loads(cache_manifest.read_text())
adata = sc.read_h5ad(auc_checkpoint)
# Preserve original plotting metadata (including annotation palettes) without
# loading the full expression matrix, raw counts, layers or embeddings.
with h5py.File(full_checkpoint, 'r') as checkpoint:
    adata.uns = read_elem(checkpoint['uns'])
auc_mtx = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
auc_mtx = auc_mtx.loc[record["auc_index"]]
auc_mtx.index.name = record["auc_index_name"]
auc_mtx.columns.name = record["auc_columns_name"]
regulons = read_regulons()
# Preserve the original RNG continuation for plots that use random jitter.
rng = record["numpy_random_state"]
np.random.set_state((rng[0], np.asarray(rng[1], dtype=np.uint32), rng[2], rng[3], rng[4]))
def nested_tuple(value):
    return tuple(nested_tuple(item) for item in value) if isinstance(value, list) else value
random.setstate(nested_tuple(record["python_random_state"]))

regulon_cols = [c for c in adata.obs.columns if c.startswith("Regulon(")]
plot_regulon_cols = regulon_cols
if len(regulon_cols) > 0 and top_regulons > 0 and len(regulon_cols) > top_regulons:
  regulon_activity_var = adata.obs[regulon_cols].var().sort_values(ascending=False)
  plot_regulon_cols = regulon_activity_var.head(top_regulons).index.tolist()
logger.info(
    f"Using {len(plot_regulon_cols)} regulon(s) for dotplot/violin visualization "
    f"(pyscenic_top_regulons={top_regulons if top_regulons > 0 else 'all'})"
)
celltype_col = args.celltype if args.celltype in adata.obs else None
if celltype_col is None:
  cat_cols = [c for c in adata.obs.columns if str(adata.obs[c].dtype) == "category"]
  celltype_col = cat_cols[0] if len(cat_cols) > 0 else args.celltype

if len(regulon_cols) > 0:
  log_step(logger, 3, 5, "saving regulon activity plots")
  for path, targets in figure_pages(f"{args.sample_id}_dotplot_regulons.png", plot_regulon_cols, 20):
    sc.pl.dotplot(adata, targets, groupby=celltype_col, standard_scale='var', dot_min=0.02,
                  figsize=(14, max(4, min(16, adata.obs[celltype_col].nunique() * 0.4 + 2))), show=False)
    plt.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

if len(regulon_cols) > 0:
  violin_targets = plot_regulon_cols if top_regulons > 0 else regulon_cols
  for path, targets in figure_pages(f"{args.sample_id}_violin_regulons.png", violin_targets, 20):
    ncols = min(4, len(targets))
    nrows = (len(targets) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), squeeze=False)
    for ax, key in zip(axes.flat, targets):
      sc.pl.violin(adata, keys=key, groupby=celltype_col, rotation=90, ax=ax, show=False)
      ax.set_title(key, fontsize=9)
      ax.set_ylabel('AUCell')
      ax.set_xlabel('')
      ax.tick_params(axis='x', labelsize=7)
    for ax in list(axes.flat)[len(targets):]:
      ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

plot_auc_heatmap_scanpy(auc_mtx, adata, groupby=celltype_col, top_genes=30, outputs=os.path.join(output_dir, f"{args.sample_id}_auc_heatmap.png"))

log_step(logger, 4, 5, "saving AUCell and regulon result tables")
auc_mtx.to_csv(os.path.join(output_dir,f'{args.sample_id}.auc.csv'))

# Save Regulon-Gene relationships for reproducibility
regulon_list = []
for reg in regulons:
    for gene in reg.genes:
        regulon_list.append({"Regulon": reg.name, "Gene": gene})
pd.DataFrame(regulon_list).to_csv(os.path.join(output_dir, f"{args.sample_id}_regulon_genes.csv"), index=False)


df_auc = auc_mtx.copy()
df_auc[celltype_col] = adata.obs.loc[df_auc.index, celltype_col].values
df_auc.groupby(celltype_col).mean().to_csv(os.path.join(output_dir,f"{args.sample_id}_auc_mean_by_celltype.csv"))

if celltype_col in adata.obs:
  rss = regulon_specificity_scores(auc_mtx, adata.obs[celltype_col])
  rss.to_csv(os.path.join(output_dir,f"{args.sample_id}_rss.csv"))
  celltypes = adata.obs[celltype_col].value_counts().index.tolist()
  if len(celltypes) > 0:
    # Save top regulons for each cell type for reproducibility
    top_regs_list = []
    for ct in celltypes:
        top_n = rss.loc[ct].sort_values(ascending=False).head(10)
        for rank, (reg, score) in enumerate(top_n.items(), 1):
            top_regs_list.append({"CellType": ct, "Rank": rank, "Regulon": reg, "RSS": score})
    pd.DataFrame(top_regs_list).to_csv(os.path.join(output_dir, f"{args.sample_id}_rss_top10.csv"), index=False)

    # Dynamic layout calculation
    ncols = 5  # Fixed at 5 columns for better readability if many types
    nrows = math.ceil(len(celltypes) / ncols)

    # Adjust figure size dynamically
    fig_width = 4 * ncols
    fig_height = 3.5 * nrows

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), dpi=120)
    axes = np.array(axes).reshape(-1)

    for i, ct in enumerate(celltypes):
      # Skip if data is invalid (e.g. all NaNs)
      if rss.loc[ct].isna().all():
        logger.warning(f"RSS data for {ct} contains only NaNs; skip plot")
        continue

      # Fill NaNs with 0 to prevent plotting errors
      rss_filled = rss.fillna(0)
      plot_rss(rss_filled, ct, ax=axes[i])

      axes[i].set_xlabel('')
      axes[i].set_ylabel('')
      # Add title to identify cell type clearly
      axes[i].set_title(ct, fontsize=10)

    for j in range(len(celltypes), len(axes)):
      axes[j].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,f"{args.sample_id}_rss.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

df_obs = adata.obs
if len(regulon_cols) > 0:
  df_scores = df_obs[regulon_cols + [celltype_col]].copy()
  std = df_obs[regulon_cols].std()
  valid_cols = std[std > 0].index.tolist()
  if len(valid_cols) > 0:
    df_scores = df_scores[valid_cols + [celltype_col]]
    df_z = (df_scores.groupby(celltype_col).mean() - df_obs[valid_cols].mean()) / df_obs[valid_cols].std()

    # Save full Z-score matrix for reproduction
    df_z.to_csv(os.path.join(output_dir, f"{args.sample_id}_zscore_matrix.csv"))

    # Plot ALL regulons sorted by mean absolute Z-score activity
    top_regs = df_z.abs().mean().sort_values(ascending=False).index
    df_heatmap = df_z[top_regs]

    vmax = max(float(np.nanmax(np.abs(df_heatmap.to_numpy()))), 1e-6)
    for path, targets in figure_pages(f"{args.sample_id}_zscore_heatmap.png", list(top_regs), 40):
      fig, ax1 = plt.subplots(figsize=(max(10, len(targets) * 0.35), max(4, min(16, len(df_heatmap) * 0.5))))
      sns.heatmap(df_heatmap[targets], ax=ax1, cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
                  linewidths=0, cbar_kws={'label': 'Z-score'})
      ax1.set_ylabel('')
      ax1.set_xlabel('')
      fig.savefig(path, dpi=300, bbox_inches='tight')
      plt.close(fig)
  violin_targets = []
  if "df_z" in locals():
    if top_regulons > 0:
      violin_targets = df_z.abs().mean().sort_values(ascending=False).head(top_regulons).index.tolist()
    else:
      violin_targets = df_z.abs().mean().sort_values(ascending=False).index.tolist()
  else:
    violin_targets = plot_regulon_cols if top_regulons > 0 else regulon_cols
  if len(violin_targets) > 0:
    medians = adata.obs.groupby(celltype_col, observed=True)[violin_targets].median().to_numpy()
    for path, targets in figure_pages(f"{args.sample_id}_stacked_violin.png", violin_targets, 20):
      fig = sc.pl.stacked_violin(adata, targets, groupby=celltype_col, return_fig=True, show=False,
                                 figsize=(max(8, len(targets) * 0.5 + 3), max(4, min(16, adata.obs[celltype_col].nunique() * 0.4 + 2))),
                                 vmin=float(np.nanmin(medians)), vmax=float(np.nanmax(medians)))
      fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
      plt.close('all')

logger.info("pySCENIC plotting completed using saved AUCell scores")
