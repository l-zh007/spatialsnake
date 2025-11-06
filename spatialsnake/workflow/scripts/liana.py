import os
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
import spatialdata_plot as splt
import pandas as pd
import scanpy as sc
import decoupler as dc
import liana as li
from matplotlib import pyplot as plt
plt.rcParams['figure.dpi'] = 50
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean
import argparse
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--bandwidth', type=float, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--cutoff', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--resource_name', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--celltype', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--expr_prop', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sp_cellchat', type=str, required=False,
                   help='Path for the output zarr file')                   
args = parser.parse_args()                   
           

output_dir=os.path.dirname(args.output_zarr_path)

def spatial_ccc_analysis(output_dir):
  plot, _ = li.ut.query_bandwidth(coordinates=adata.obsm['spatial'], start=0, end=500, interval_n=20)
  plot.save(os.path.join(output_dir,'bandwidth_plot.png'), width=8, height=6, dpi=300)
  li.ut.spatial_neighbors(adata, bandwidth=args.bandwidth, cutoff=args.cutoff, kernel='gaussian', set_diag=True)
  plot2=li.pl.connectivity(adata, idx=200, size=0.05, figure_size=(6, 5))
  plot2.save(os.path.join(output_dir,'connectivity.png'), width=8, height=6, dpi=300)

  lrdata = li.mt.bivariate(adata,
                resource_name=args.resource_name, # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                nz_prop=0.05, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True)
  
  fig, axes = plt.subplots(2, 2, figsize=(12, 10))
  axes = axes.flatten()
  top_interactions_list = [lrdata.var.sort_values("morans", ascending=False).index[0],lrdata.var.sort_values("mean", ascending=False).index[0], lrdata.var.sort_values("std", ascending=False).index[0],lrdata.var.sort_values("morans_pvals", ascending=True).index[0]]
  titles = ["Highest Spatial\n(Moran's I)", "Highest Mean", "Highest Variability\n(Std)", "Most Significant"]
  for idx, (interaction, title) in enumerate(zip(top_interactions_list, titles)):
    sc.pl.spatial(lrdata, color=interaction, size=1.4, vmax=1, cmap='magma', ax=axes[idx], show=False, title=f'{title}\n{interaction}')
  plt.tight_layout()
  plt.savefig("top_interactions_simple.png", dpi=300, bbox_inches='tight')
  plt.close()
  
  
  li.multi.nmf(lrdata, n_components=None, inplace=True, random_state=0, max_iter=200, verbose=True)
  lr_loadings = li.ut.get_variable_loadings(lrdata, varm_key='NMF_H').set_index('index')
  factor_scores = li.ut.get_factor_scores(lrdata, obsm_key='NMF_W')
  nmf = sc.AnnData(X=lrdata.obsm['NMF_W'],
                 obs=lrdata.obs,
                 var=pd.DataFrame(index=lr_loadings.columns),
                 uns=lrdata.uns,
                 obsm=lrdata.obsm)
  sc.pl.spatial(nmf, color=[*nmf.var.index, None], size=1.4, ncols=2)
  plt.savefig(os.path.join(output_dir,"factor_scores.png"),dpi=300,bbox_inches='tight')
  plt.close()
  
  lrdata.write(os.path.join(output_dir,args.output_zarr_path))
  nmf.write(os.path.join(output_dir,'nmf_spatial_factors.h5ad'))
  lrdata.var.to_csv(os.path.join(output_dir,'interaction_statistics.csv'))
  lr_loadings.to_csv(os.path.join(output_dir,'nmf_factor_loadings.csv'))
  factor_scores.to_csv(os.path.join(output_dir,'nmf_factor_scores.csv'))

def stable_ccc_analysis(output_dir):
  cellphonedb(adata,
            groupby=args.celltype, 
            resource_name=args.resource_name,
            expr_prop=args.expr_prop, #0.1
            verbose=True, key_added='cpdb_res')
  if not source_labels and target_labels:
    source_labels = adata.uns['cpdb_res']['source'].unique().tolist()
    target_labels = adata.uns['cpdb_res']['target'].unique().tolist()
  if not int_celltype:
    all_celltypes = source_labels[0]
  dotplot_fig=li.pl.dotplot(adata = adata, 
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=source_labels[:2],
              target_labels=target_labels[-2:],
              figure_size=(12, 10),
              # finally, since cpdbv2 suggests using a filter to FPs
              # we filter the pvals column to <= 0.05
              filter_fun=lambda x: x['cellphone_pvals'] <= 0.05,
              uns_key='cpdb_res')
  dotplot_fig.save(os.path.join(output_dir,'dotplot.png'),dpi=300)
  my_plot = li.pl.tileplot(adata = adata,
                         fill='means',
                         label='props',
                         label_fun=lambda x: f'{x:.2f}',
                         top_n=10, 
                         orderby='cellphone_pvals',
                         orderby_ascending=True,
                         source_labels=source_labels[:2],
                         target_labels=target_labels[-2:],
                         uns_key='cpdb_res', # NOTE: default is 'liana_res'
                         source_title='Ligand',
                         target_title='Receptor',
                         figure_size=(10, 8))
  my_plot.save(os.path.join(output_dir,'tileplot.png'),dpi=300)
  circle=li.pl.circle_plot(adata,
                  groupby=args.celltype,
                  score_key='lr_means',
                  inverse_score=True,
                  source_labels=all_celltypes[0],
                  filter_fun=lambda x: x['cellphone_pvals'] <= 0.05,
                  pivot_mode='counts', # NOTE: this will simply count the interactions, 'mean' is also available
                  figure_size=(10, 10),
                  uns_key='cpdb_res')
  plt.savefig(os.path.join(output_dir,'circle.png'), dpi=300, bbox_inches='tight', facecolor='white')
  plt.close()
  adata.write(os.path.join(output_dir,args.output_zarr_path))    
if os.path.splitext(args.input_dir)[1].lower()==".h5ad":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
if args.sp_cellchat=="False":
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  systems=[]
  if len(valid_coord_systems)>1:
    for i in range(len(valid_coord_systems)):
        if args.coord_type in valid_coord_systems[i]:
          systems.append(valid_coord_systems[i])
  else:
    systems=valid_coord_systems
  adata = to_legacy_anndata(concatenated_sdata, include_images=True, coordinate_system=systems[0])
  adata.raw = adata.copy()
  cord=pd.DataFrame(data=adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
  adata.obs['array_row']=cord["x"]
  adata.obs['array_col']=cord["y"]
  spatial_ccc_analysis(output_dir)
else:
  stable_ccc_analysis(output_dir)
    
                   
                   
                   
                   
                   
                   
                   
                   
