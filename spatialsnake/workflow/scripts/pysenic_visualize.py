import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")
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
from pyscenic.aucell import aucell
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from spatialsnake.workflow.function.plot import plot_auc_heatmap_scanpy
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step
import argparse
logger = setup_logger("pyscenic_visualize")
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
args = parser.parse_args()
output_dir = os.path.dirname(args.outputs)
os.makedirs(output_dir,exist_ok=True)
top_regulons = max(0, int(args.top_regulons))
min_regulon_genes = max(1, int(args.min_regulon_genes))

log_step(logger, 1, 5, "loading expression data and regulons")
if os.path.splitext(args.input_dir)[1].lower()==".h5ad":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]

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

log_step(logger, 2, 5, "calculating AUCell matrix")
exp_matrix = adata.X
if issparse(exp_matrix):
  nnz = exp_matrix.nnz
  dense_bytes = adata.n_obs * adata.n_vars * np.dtype(np.float32).itemsize
  logger.info(
      f"Sparse expression matrix detected with shape={exp_matrix.shape}, nnz={nnz}, "
      f"estimated dense float32 memory={dense_bytes / 1024**3:.2f} GB"
  )
  if dense_bytes > 8 * 1024**3:
    raise MemoryError(
        "pySCENIC visualization requires dense AUCell ranking input and the current matrix is estimated to exceed 8 GB "
        "when densified. Please reduce object size before visualization."
    )
  exp_matrix = exp_matrix.toarray()
else:
  logger.info(f"Dense expression matrix detected with shape={exp_matrix.shape}")
gene_names = adata.var_names
cell_names = adata.obs_names
exp_df = pd.DataFrame(exp_matrix,index=cell_names,columns=gene_names)
auc_mtx = aucell(exp_df, regulons, num_workers=args.num_workers)

add_scenic_metadata(adata, auc_mtx, regulons)
adata.write(os.path.join(output_dir,f"{args.sample_id}.h5ad"))

aucell_adata = sc.AnnData(X=auc_mtx.loc[adata.obs_names].values)
aucell_adata.obs = adata.obs
aucell_adata.var_names = auc_mtx.columns
aucell_adata.write(os.path.join(output_dir,f"{args.sample_id}_aucell.h5ad"))

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
  sc.pl.dotplot(adata, plot_regulon_cols, groupby=celltype_col, standard_scale='var', dot_min=0.02, show=False)
  plt.savefig(os.path.join(output_dir,f"{args.sample_id}_dotplot_regulons.png"), dpi=300, bbox_inches='tight', facecolor='white')
  plt.close()

if len(regulon_cols) > 0:
  violin_targets = plot_regulon_cols if top_regulons > 0 else regulon_cols
  sc.pl.violin(adata, keys=violin_targets, groupby=celltype_col, multi_panel=True, show=False)
  plt.savefig(os.path.join(output_dir,f"{args.sample_id}_violin_regulons.png"), dpi=300, bbox_inches='tight', facecolor='white')
  plt.close()

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
    
    # Adjust figure size based on number of regulons
    n_cols = df_heatmap.shape[1]
    n_rows = df_heatmap.shape[0]
    fig_width = max(12, n_cols * 0.25)
    fig_height = max(8, n_rows * 0.5)
    
    fig, ax1 = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    # Use linewidths=0 for cleaner look with many columns
    sns.heatmap(df_heatmap, ax=ax1, cmap="vlag", center=0, linewidths=0, cbar_kws={'label': 'Z-score'})
    ax1.set_ylabel('')
    ax1.set_xlabel('')
    fig.savefig(os.path.join(output_dir,f"{args.sample_id}_zscore_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.close()
  violin_targets = []
  if "df_z" in locals():
    if top_regulons > 0:
      violin_targets = df_z.abs().mean().sort_values(ascending=False).head(top_regulons).index.tolist()
    else:
      violin_targets = df_z.abs().mean().sort_values(ascending=False).index.tolist()
  else:
    violin_targets = plot_regulon_cols if top_regulons > 0 else regulon_cols
  if len(violin_targets) > 0:
    fig = sc.pl.stacked_violin(adata, violin_targets, groupby=celltype_col, return_fig=True, show=False)
    fig.savefig(os.path.join(output_dir, f"{args.sample_id}_stacked_violin.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
log_step(logger, 5, 5, f"writing pySCENIC loom output to {args.outputs}")
adata.write_loom(args.outputs)
logger.info("pySCENIC visualization module completed")
