import os
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import gc
import argparse
from spatialdata.datasets import blobs_annotating_element
from spatialdata.transformations import Affine, set_transformation
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
import spatialdata_plot as splt
import pandas as pd
import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
from cytoolz import compose
import operator as op
import seaborn as sns
import anndata as ad
from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
import os, time, random, gc
import anndata
from anndata import AnnData
import numpy as np
import pandas as pd
import warnings
from banksy_utils.color_lists import spagcn_color
warnings.filterwarnings("ignore")

concatenated_sdata = spd.read_zarr("./results/useful_results/Lesional_1.zarr")
print(concatenated_sdata)
for table in concatenated_sdata.tables.keys():
  table=table
  adata = concatenated_sdata[table]
  
print(adata)

coord_keys = ('array_row','array_col','coord_xy')
adata.obsm[coord_keys[2]] = np.vstack((adata.obs[coord_keys[0]].values,adata.obs[coord_keys[1]].values)).T

print(adata.obs)

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10,10)
cmap = plt.get_cmap('tab20')
sc.pl.scatter(adata, x='array_row', y='array_col', color='celltype',  color_map=cmap, title=f"Tissue", size = 5)#show = False, size = 5)
plt.show()
plt.savefig('fig/scatter.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()




from banksy.main import median_dist_to_nearest_neighbour
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.main import concatenate_all

k_geom = 15  # only for fixed type
max_m = 1  # azumithal transform up to kth order
nbr_weight_decay = "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"
resolutions = None  # clustering resolution for leiden algorithm
max_labels = 8 # Number of clusters for tissue segmentation
pca_dims = [20]  # Dimensionality in which PCA reduces to
lambda_list = [0.8]
file_path = "fig"


nbrs = median_dist_to_nearest_neighbour(adata, key=coord_keys[2])
banksy_dict = initialize_banksy(adata,
                                coord_keys,
                                k_geom,
                                nbr_weight_decay=nbr_weight_decay,
                                max_m=max_m,
                                plt_edge_hist=False,
                                plt_nbr_weights=True,
                                plt_agf_angles=False,
                                plt_theta=False
                                )

banksy_dict, banksy_matrix = generate_banksy_matrix(adata,
                                                    banksy_dict,
                                                    lambda_list,
                                                    max_m)


from banksy_utils.umap_pca import pca_umap

pca_umap(banksy_dict,
         pca_dims = pca_dims,
         add_umap = True
         )




from banksy.cluster_methods import run_Leiden_partition

banksy_df, max_num_labels = run_Leiden_partition(
    banksy_dict,
    resolutions,
    num_nn = 50,
    num_iterations = -1,
    partition_seed = 1234,
    match_labels = True,
    max_labels = max_labels,
)

from banksy.plot_banksy import plot_results

c_map =  'tab20' # specify color map
weights_graph =  banksy_dict['scaled_gaussian']['weights'][1]

plot_results(
    banksy_df,
    weights_graph,
    c_map,
    match_labels = True,
    coord_keys = coord_keys,
    max_num_labels  =  max_num_labels, 
    save_path = os.path.join(file_path, 'BANKSY-Results'),
    save_fig = False, # Save Spatial Plot Only
    save_fullfig = True, # Save Full Plot
    dataset_name = f"CODEX",
    save_labels=True
)
banksy_df.to_csv(os.path.join(file_path, f"BANKSY.csv"))


# Add nonspatial clustering
nonspatial_dict = {"nonspatial" : {0.0: {"adata": concatenate_all([adata.X], 0, adata=adata), } } }

pca_umap(nonspatial_dict, pca_dims = pca_dims, add_umap = True )

from banksy.cluster_methods import run_Leiden_partition

nonspatial_df, max_num_labels = run_Leiden_partition(
    nonspatial_dict,
    resolutions,
    num_nn = 50,
    num_iterations = -1,
    partition_seed = 1234,
    match_labels = True,
    max_labels = max_labels,
)

from banksy.plot_banksy import plot_results

c_map =  'tab20' # specify color map

plot_results(
    nonspatial_df,
    weights_graph,
    c_map,
    match_labels = True,
    coord_keys = coord_keys,
    max_num_labels  =  max_num_labels, 
    save_path = os.path.join(file_path, 'BANKSY-Results'),
    save_fig = False, # Save Spatial Plot Only
    save_fullfig = True, # Save Full Plot
    dataset_name = f"CODEX-",
    save_labels=True
)


from sklearn.metrics import adjusted_rand_score as ari, adjusted_mutual_info_score as ami
from sklearn.metrics import matthews_corrcoef as mcc
# See the visualize the communities that we want to detect
adata.obs['celltype']


banksy_spatial_clusters = banksy_df.labels[banksy_df.index[0]]
banksy_spatial_clusters.dense

nonspatial_clusters = nonspatial_df.labels[nonspatial_df.index[0]]
nonspatial_clusters.dense


def calculate_metrics(cluster_labels, annotated_labels):
    # A custom function to calcualte all metrics
    ari_score  = ari(cluster_labels, annotated_labels)
    ami_score =   ami(cluster_labels, annotated_labels)

    if isinstance(annotated_labels.dtype, pd.CategoricalDtype):
        print("Converting annotations to required 'int' type for computing MCC")
        annotated_labels = annotated_labels.cat.codes

    mcc_score =  mcc(cluster_labels,annotated_labels )
    return ari_score, ami_score, mcc_score

nonspatial_ari, nonspatial_ami, nonspatial_mcc = calculate_metrics(nonspatial_clusters.dense, adata.obs['celltype'])


banksy_ari, banksy_ami, banksy_mcc = calculate_metrics(banksy_spatial_clusters.dense, adata.obs['celltype'])


def bar_plot(metrics, methods):
    ''' Custom function to generate bar chart comparing metrices of labels produced by different methods'''
    fig, ax = plt.subplots(figsize=(8,8),layout='constrained')
    x = np.arange(len(methods))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0
    # method is banksy, nonspatial,
    # metric is ari, ami, mcc
    for method, metric  in metrics.items():
        offset = width * multiplier
        print(metric)
        rects = ax.bar(x + offset, metric, width, label=method)
        ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Metrices', fontsize=18)
    ax.set_title("Similarity between BANKSY labels and annotated communities", fontsize=20)
    ax.set_xticks(x + width, methods, fontsize=16)
    ax.legend(loc='upper left', fontsize=18)
    fig.show()
    plt.savefig(os.path.join(file_path, 'BANKSY-Results/bar.png'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

### Plot the similarity between BANKSY labels and annotated communities
methods = ('Non-spatial labels', 'BANKSY labels')
metrics = {
    'Adjusted Rand Index' : (nonspatial_ari, banksy_ari),
    'Adjusted Mutual Information': (nonspatial_ami, banksy_ami),
    'Matthew Correlation Coefficient': (nonspatial_mcc,banksy_mcc ),
}
bar_plot(metrics, methods)

