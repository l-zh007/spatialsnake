import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")
import spatialdata as spd
import scanpy as sc
import spatialdata_plot as splt
import pandas as pd
import re, pickle, math
from functools import partial
from collections import OrderedDict
from cytoolz import compose
import anndata as ad
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
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from spatialsnake.workflow.function.plot import plot_auc_heatmap_scanpy
import argparse
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
args = parser.parse_args()
output_dir = os.path.dirname(args.outputs)
os.makedirs(output_dir,exist_ok=True)

if args.types=="slide_seq" or os.path.splitext(args.input_dir)[1].lower()==".h5ad":
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
  regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0)])))
  regulons = list(map(lambda r: r.rename(r.transcription_factor), regulons))

exp_matrix = adata.X
if issparse(exp_matrix):
  exp_matrix = exp_matrix.toarray()
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
celltype_col = args.celltype if args.celltype in adata.obs else None
if celltype_col is None:
  cat_cols = [c for c in adata.obs.columns if str(adata.obs[c].dtype) == "category"]
  celltype_col = cat_cols[0] if len(cat_cols) > 0 else args.celltype

if len(regulon_cols) > 0:
  sc.pl.dotplot(adata, regulon_cols, groupby=celltype_col, standard_scale='var', dot_min=0.02, show=False)
  plt.savefig(os.path.join(output_dir,f"{args.sample_id}_dotplot_regulons.png"), dpi=300, bbox_inches='tight', facecolor='white')
  plt.close()

if len(regulon_cols) > 0:
  n_show = min(12, len(regulon_cols))
  sc.pl.violin(adata, keys=regulon_cols[:n_show], groupby=celltype_col, multi_panel=True, show=False)
  plt.savefig(os.path.join(output_dir,f"{args.sample_id}_violin_regulons.png"), dpi=300, bbox_inches='tight', facecolor='white')
  plt.close()

plot_auc_heatmap_scanpy(auc_mtx, adata, groupby=celltype_col, top_genes=20, outputs=os.path.join(output_dir, f"{args.sample_id}_auc_heatmap.png"))

auc_mtx.to_csv(os.path.join(output_dir,f'{args.sample_id}.auc.csv'))

df_auc = auc_mtx.copy()
df_auc[celltype_col] = adata.obs.loc[df_auc.index, celltype_col].values
df_auc.groupby(celltype_col).mean().to_csv(os.path.join(output_dir,f"{args.sample_id}_auc_mean_by_celltype.csv"))

if celltype_col in adata.obs:
  rss = regulon_specificity_scores(auc_mtx, adata.obs[celltype_col])
  rss.to_csv(os.path.join(output_dir,f"{args.sample_id}_rss.csv"))
  celltypes = adata.obs[celltype_col].value_counts().index.tolist()
  celltypes = celltypes[:8]
  if len(celltypes) > 0:
    ncols = 4 if len(celltypes) > 4 else len(celltypes)
    nrows = math.ceil(len(celltypes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), dpi=120)
    axes = np.array(axes).reshape(-1)
    for i, ct in enumerate(celltypes):
      plot_rss(rss, ct, ax=axes[i])
      axes[i].set_xlabel('')
      axes[i].set_ylabel('')
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
    top_regs = df_z.abs().mean().sort_values(ascending=False).head(30).index
    df_heatmap = df_z[top_regs]
    fig, ax1 = plt.subplots(1, 1, figsize=(12, 8))
    sns.heatmap(df_heatmap, ax=ax1, cmap="vlag", center=0, linewidths=.6, linecolor='white')
    ax1.set_ylabel('')
    fig.savefig(os.path.join(output_dir,f"{args.sample_id}_zscore_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.close()
  violin_targets = []
  if "df_z" in locals():
    violin_targets = df_z.abs().mean().sort_values(ascending=False).head(20).index.tolist()
  else:
    violin_targets = regulon_cols[:20]
  if len(violin_targets) > 0:
    fig = sc.pl.stacked_violin(adata, violin_targets, groupby=celltype_col, return_fig=True, show=False)
    fig.savefig(os.path.join(output_dir, f"{args.sample_id}_stacked_violin.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
adata.write_loom(args.outputs)


