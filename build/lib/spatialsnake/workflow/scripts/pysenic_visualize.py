import spatialdata as spd
import scanpy as sc
import spatialdata_plot as splt
import pandas as pd
import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
from cytoolz import compose
import anndata as ad
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
import numpy as np
from matplotlib import pyplot as plt
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


if types=="slide_seq" or os.path.splitext(args.input_dir)[1].lower()==".h5ad":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]

def derive_regulons(code, folder):
    motifs = load_motifs(args.regulons)
    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0)])))
    regulons = list(map(lambda r: r.rename(r.transcription_factor), regulons))
    with open(os.path.join(folder, '{}.regulons.dat'.format(code)), 'wb') as f:
        pickle.dump(regulons, f)
    return os.path.join(folder, '{}.regulons.dat'.format(code))


paths=derive_regulons(args.sample_id, output_dir)
with open(paths, 'rb') as f:
    regulons = pickle.load(f)

exp_matrix = adata.X
if issparse(exp_matrix):
  exp_matrix = exp_matrix.toarray()
gene_names = adata.var_names
cell_names = adata.obs_names
exp_df = pd.DataFrame(exp_matrix,index=cell_names,columns=gene_names)
auc_mtx = aucell(exp_df, regulons, num_workers=args.num_workers)
add_scenic_metadata(adata, auc_mtx, regulons)
sc.tl.umap(adata, use_rep='X_aucell')
aucell_adata = sc.AnnData(X=auc_mtx.sort_index())
aucell_adata.obs = adata.obs


gene_names = []
for col in adata.obs.columns:
  match = re.search(r'Regulon\((.*?)\)', col)
  if match:
    gene_names.append(match.group(1))
    

sc.pl.umap(adata, color=gene_names,title=gene_names, ncols=2, use_raw=False)
plt.savefig('senic_tsne.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()


sc.pl.dotplot(adata, gene_names, groupby=args.celltype)
plt.savefig('senic_dotplot.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

plot_auc_heatmap_scanpy(auc_mtx, adata, groupby=args.celltype, top_genes=20,outputs=args.outputs)

auc_mtx.to_csv(os.path.join(output_dir,f'{args.sample_id}.auc.csv'))
concatenated_sdata[table]=adata
concatenated_sdata.write(os.path.join(output_dir,"scenic_results.zarr"))














