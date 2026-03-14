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

#adata = sc.read_loom(os.path.join(output_dir,"breast_cancer2.loom"))
#print(adata.obs_names)
#import pandas as pd
#df = pd.read_csv("../spatial_project/FU.csv", index_col=0)
#adata.obs_names = df.index.astype(str)
#adata.var_names = df.columns.astype(str)
#adata.X = df.values
#print(adata.obs_names)
#df = pd.read_csv("../spatial_project/celltype_FU.csv")
#print(df)
#df = df.set_index("cell_id")
#adata.obs["celltype"] = df.loc[adata.obs_names, "celltype"].values
#print(adata)



if load_regulons is not None:
  regulons = load_regulons("results/pysenic_results/FU_pysenic/breast_cancer2.regulons.csv")
else:
  motifs = load_motifs("results/pysenic_results/FU_pysenic/breast_cancer2.regulons.csv")
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

plot_auc_heatmap_scanpy(auc_mtx, adata, groupby=celltype_col, top_genes=30, outputs=os.path.join(output_dir, f"{args.sample_id}_auc_heatmap.png"))

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
        print(f"Warning: RSS data for {ct} contains only NaNs. Skipping plot.")
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
    violin_targets = df_z.abs().mean().sort_values(ascending=False).head(20).index.tolist()
  else:
    violin_targets = regulon_cols[:20]
  if len(violin_targets) > 0:
    fig = sc.pl.stacked_violin(adata, violin_targets, groupby=celltype_col, return_fig=True, show=False)
    fig.savefig(os.path.join(output_dir, f"{args.sample_id}_stacked_violin.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
adata.write_loom(args.outputs)
