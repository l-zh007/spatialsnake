import os
import numpy as np
import pandas as pd
import scanpy as sc
from spatialsnake.workflow.banksy.main import median_dist_to_nearest_neighbour
from spatialsnake.workflow.banksy.initialize_banksy import initialize_banksy
from spatialsnake.workflow.banksy.embed_banksy import generate_banksy_matrix
from spatialsnake.workflow.banksy.main import concatenate_all
from spatialsnake.workflow.banksy_utils.umap_pca import pca_umap
from spatialsnake.workflow.banksy.cluster_methods import run_Leiden_partition
from spatialsnake.workflow.banksy.plot_banksy import plot_results
import matplotlib.pyplot as plt
import spatialdata as spd
import scanpy.external as sce
import matplotlib.pyplot as plt
from spatialdata.datasets import blobs_annotating_element
from spatialdata.transformations import Affine, set_transformation
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
import spatialdata_plot as splt
import argparse
############################clustering################################
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--tsene', type=bool, required=False,
                   help='Path for the output zarr file') 
parser.add_argument('--MIN_DIST', type=float, required=False,
                   help='Path for the output zarr file')                    
parser.add_argument('--SPREAD', type=float, required=False,
                   help='Path for the output zarr file') 

parser.add_argument('--RES', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--cluster_algorithm', type=str, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--n_clusters', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--k_geom', type=float, required=False,
                   help='Path for the output zarr file')                   
parser.add_argument('--max_m', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--lambda_list', type=float, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--nbr_weight_decay', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--n_comps', type=int, required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()

type=args.type

dir_path=os.path.dirname(args.output_zarr_path)

k_geom = args.k_geom
max_m = args.max_m
nbr_weight_decay = args.nbr_weight_decay
resolutions = [args.RES]
pca_dims = [args.n_comps]
lambda_list = [args.lambda_list]

if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
  
sc.tl.umap(adata,min_dist=args.MIN_DIST, spread=args.SPREAD, random_state=0)
if args.tsene==True:
  sc.tl.tsne(adata, n_pcs=50, perplexity=30)
sdata = to_legacy_anndata(concatenated_sdata, include_images=True, coordinate_system="downscale_to_hires")
cord=pd.DataFrame(data=sdata.obsm['spatial'],index=sdata.obs_names,columns=['x','y'])
sdata.obs['array_row']=cord["x"]
sdata.obs['array_col']=cord["y"]
nbrs = median_dist_to_nearest_neighbour(sdata, key="spatial")
coord_keys = ('array_row', 'array_col', 'spatial')



banksy_dict = initialize_banksy(
    sdata,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=True,
    plt_nbr_weights=True,
    plt_agf_angles=False,
    plt_theta=True)

banksy_dict, banksy_matrix = generate_banksy_matrix(
                                sdata,
                                banksy_dict,
                                lambda_list,
                                max_m)
pca_umap(
    banksy_dict,
    pca_dims = pca_dims,
    add_umap = True,
    plt_remaining_var = False,
)
results_df, max_num_labels = run_Leiden_partition(
    banksy_dict,
    resolutions,
    num_nn = 50,
    num_iterations = -1,
    partition_seed = 12345,
    match_labels = True)
adata.obs["clusters"]=banksy_matrix.obs["labels_scaled_gaussian_pc20_nc0.20_r0.50"]
adata.obs['clusters'] = adata.obs['clusters'].astype('category')




sc.pl.umap(adata,color=[
    "total_counts",
    "n_genes_by_counts",
    "clusters"],wspace=0.4)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}UMAP.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()
    

if args.tsene==True:
  sc.pl.tsne(adata, color=[
      "total_counts",
      "n_genes_by_counts",
      "clusters"], ncols=3, show=True)
  plt.savefig(
      os.path.join(dir_path,f"{args.sample_id}tsene.png"),
      dpi=300,
      bbox_inches='tight')
  plt.show()
  plt.close()



    
sample_names = adata.obs["region"].unique()
plt.imshow(pd.crosstab(adata.obs["region"], adata.obs["clusters"]), cmap='hot', interpolation='nearest')
plt.title("Cell Distribution Across Clusters")
plt.xlabel("Cluster")
plt.yticks(range(len(sample_names)), sample_names)
plt.savefig(
      os.path.join(dir_path, f"{args.sample_id}Cell_Distribution_Across_Clusters.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()




if type!="slide_seq":
    # if type=="visium":
    #   adata.obs['sopt_id'] = adata.obs['spot_id'].astype(str)
    #   adata.obs['region'] = adata.obs['region'].astype('category')
    # else:
    adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
    adata.obs['region'] = adata.obs['region'].astype('category')
    concatenated_sdata[table]=adata
    concatenated_sdata.write(
      os.path.join(args.output_zarr_path),
      overwrite=True)
else:  
    adata.write(os.path.join(args.output_zarr_path))












