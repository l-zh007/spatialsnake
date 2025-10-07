import os
# theano
# os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
# def set_theano_device():
#     try:
#         import theano
#         theano.config.device = 'cuda'
#         theano.config.floatX = 'float32'
#         return 'cuda'
#     except (ImportError, theano.gpuarray.GpuArrayException):
#         print("GPU不可用，使用CPU")
#         return 'cpu'
# device = set_theano_device()
import sys
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import spatialdata_io
import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
import argparse
import matplotlib as mpl
import cell2location
from matplotlib import rcParams

from typing import Literal
from spatialdata.datasets import blobs_annotating_element
from spatialdata.transformations import Affine, set_transformation
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
import argparse
import seaborn as sns
from scanpy.get import obs_df
from matplotlib.pyplot import get_cmap
import itertools
import logging
from matplotlib import use
from scipy.sparse import issparse
from scvi import REGISTRY_KEYS
from typing import Union, List, Optional, Iterable, Sequence, Dict
from matplotlib.axes import Axes
from anndata import AnnData
import matplotlib




















parser = argparse.ArgumentParser()
parser.add_argument("--input_spatial", required=True)
parser.add_argument("--input_singlecell", required=True)
parser.add_argument("--output_dir_zarr",  required=True)
parser.add_argument("--max_epochs_reference",type=int, required=True)
parser.add_argument("--remove_mt",type=bool, required=True)
parser.add_argument("--sample_id",type=bool, required=True)
parser.add_argument("--type",type=str, required=True)
parser.add_argument("--max_epochs_st",type=int,required=True)
parser.add_argument("--N_cells_per_location",type=int,required=True)



parser.add_argument("--labels_key_reference", default='Subset', required=False)
parser.add_argument("--batch_key_reference", default='Sample', required=False)
parser.add_argument("--cell_count_cutoff", default=15, required=False)
parser.add_argument("--cell_percentage_cutoff2", default=0.05, required=False)
parser.add_argument("--nonz_mean_cutoff", default=1.12, required=False)
parser.add_argument("--labels_key_st", default='Subset', required=False)
parser.add_argument("--batch_key_st", default='sample', required=False)
parser.add_argument("--layer_st", default=None, required=False)
parser.add_argument("--detection_alpha", default=20, required=False)
parser.add_argument("--device",required=False)
parser.add_argument("--save_models", default=True, required=False)
args = parser.parse_args()
print(args.cell_count_cutoff)
print(args.batch_key_reference)
if args.device == 'cuda':
    os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
else:
    os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32'

def cell2loc_plot_history(model,fig_path, iter_start=0, iter_end=-1, ax=None):
# Adapted from: https://github.com/BayraktarLab/cell2location/blob/master/cell2location/models/base/_pyro_mixin.py#L407
    if ax is None:
        ax = plt.gca()
    if iter_end == -1:
        iter_end = len(model.history_["elbo_train"])

    ax.plot(
        np.array(model.history_["elbo_train"].index[iter_start:iter_end]),
        np.array(model.history_["elbo_train"].values.flatten())[iter_start:iter_end],
        label="train",
    )
    ax.legend()
    ax.set_xlim(0, len(model.history_["elbo_train"]))
    ax.set_xlabel("Training epochs")
    ax.set_ylabel("-ELBO loss")
    plt.tight_layout()
    plt.savefig(fig_path)

def cell2loc_plot_QC_reference(reference_model,fig_path_reconstr, fig_path_expr,
                               summary_name: str = "means",
                               use_n_obs: int = 1000,
                               scale_average_detection: bool = True,):
# Adapted from:https://github.com/DendrouLab/panpipes/blob/main/panpipes/funcs/plotting.py

    # 1. plot reconstruction accuracy
    cell2loc_plot_QC_reconstr(model=reference_model,fig_path=fig_path_reconstr, summary_name=summary_name, use_n_obs=use_n_obs)

    # 2. plot estimated reference expression signatures (accounting for batch effect) compared to average expression in each cluster
    inf_aver = reference_model.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
    if scale_average_detection and ("detection_y_c" in list(reference_model.samples[f"post_sample_{summary_name}"].keys())):
        inf_aver = inf_aver * reference_model.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
    aver = reference_model._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
    aver = aver[reference_model.factor_names_]

    plt.hist2d(
        np.log10(aver.values.flatten() + 1),
        np.log10(inf_aver.flatten() + 1),
        bins=50,
        norm=matplotlib.colors.LogNorm(),
    )
    plt.xlabel("Mean expression for every gene in every cluster")
    plt.ylabel("Estimated expression for every gene in every cluster")
    plt.savefig(fig_path_expr)

def cell2loc_plot_QC_reconstr(model, fig_path, summary_name: str = "means", use_n_obs: int = 1000):
# Adapted from: https://github.com/BayraktarLab/cell2location/blob/master/cell2location/models/base/_pyro_mixin.py#L544

    if getattr(model, "samples", False) is False:
        raise RuntimeError("self.samples is missing, please run self.export_posterior() first")
    if use_n_obs is not None:
        ind_x = np.random.choice(
            model.adata_manager.adata.n_obs, np.min((use_n_obs, model.adata.n_obs)), replace=False
        )
    else:
        ind_x = None

    model.expected_nb_param = model.module.model.compute_expected(
        model.samples[f"post_sample_{summary_name}"], model.adata_manager, ind_x=ind_x
    )
    x_data = model.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[ind_x, :]
    if issparse(x_data):
        x_data = np.asarray(x_data.toarray())
    mu = model.expected_nb_param["mu"]
    plt.hist2d(
        np.log10(x_data.flatten() + 1),
        np.log10(mu.flatten() + 1),
        bins=50,
        norm=matplotlib.colors.LogNorm(),
    )
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("Data, log10")
    plt.ylabel("Posterior expected value, log10")
    plt.title("Reconstruction accuracy")
    plt.tight_layout()
    plt.savefig(fig_path)












output_dir=os.path.dirname(args.output_dir_zarr)

os.makedirs(os.path.join(output_dir,"figure"), exist_ok=True)

concatenated_sdata = spd.read_zarr(args.input_spatial)
for table in concatenated_sdata.tables.keys():
    table=table
    adata_vis = concatenated_sdata[table]

adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)


adata_ref = sc.read(args.input_singlecell)

adata_ref.var['SYMBOL'] = adata_ref.var.index

adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)
print(adata_ref)


if args.remove_mt:
  print("removing the MT gene")
  adata_vis.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis.var.index]
  adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]

shared_features = [f for f in adata_vis.var_names if f in adata_ref.var_names]
adata_ref = adata_ref[:, shared_features]
adata_vis = adata_vis[:, shared_features]



from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=float(args.cell_count_cutoff), cell_percentage_cutoff2=float(args.cell_percentage_cutoff2), nonz_mean_cutoff=float(args.nonz_mean_cutoff))
# filter the object
adata_ref = adata_ref[:, selected].copy()




cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        labels_key=args.labels_key_reference,  # 细胞类型标签列名
                        batch_key=args.batch_key_reference,  # 批次列名
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['Method']
                       )

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()
# 
# 
# 
# 
# 
print(args.max_epochs_reference,"###############")

mod.train(max_epochs=args.max_epochs_reference)

cell2loc_plot_history(mod, output_dir + "/figure/ELBO_sc_model.png")

adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={
        'num_samples': 1000,
        'batch_size': 2500
    }
)

if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]].copy()
else:
    inf_aver = adata_ref.var[[f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]].copy()
inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
inf_aver.to_csv(output_dir+"/Cell2Loc_inf_aver.csv")





cell2loc_plot_QC_reference(mod, output_dir + "/figure/QC_reference_reconstruction_accuracy.png", output_dir + "/figure/QC_reference_expression signatures_vs_avg_expression.png")


# if hasattr(mod, 'samples') and mod.samples is not None:
#     print("后验采样成功，开始绘制QC图")
#     mod.plot_QC()
# else:
#     print("警告：后验采样未成功，跳过QC图绘制")






# mod.train(max_epochs=250)
# 
# # mod.plot_history(20)
# # plt.savefig('./01-mod.plot_history.png')
# 
# 
# adata_ref = mod.export_posterior(
#     adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
# )

if args.save_models:
  mod.save(output_dir +"/Reference_model", overwrite=True)




# if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
#     inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
#                             for i in adata_ref.uns['mod']['factor_names']]].copy()
# else:
#     inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
#                             for i in adata_ref.uns['mod']['factor_names']]].copy()
# 
# inf_aver.columns = adata_ref.uns['mod']['factor_names']
# inf_aver.iloc[0:12, 0:12]


intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=args.batch_key_st)



model_spatial = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=args.N_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=args.detection_alpha
)
model_spatial.view_anndata_setup()


model_spatial.train(max_epochs=args.max_epochs_st,batch_size=None,train_size=1)


    
cell2loc_plot_history(model_spatial, output_dir + "/figure/ELBO_spatial_model.png")
    

adata_vis = model_spatial.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)




cell2loc_plot_QC_reconstr(model_spatial, output_dir + "/figure/QC_spatial_reconstruction_accuracy.png")

model_spatial.save(output_dir +"/Spatial_model", overwrite=True)


concatenated_sdata[table]=adata_vis
concatenated_sdata.write(output_dir +"/tem.zarr")


























# print(concatenated_sdata)
# 
# 
# adata_file = "./sp.h5ad"
# adata_vis = sc.read_h5ad(adata_file)
# mod = cell2location.models.Cell2location.load("./", adata_vis)
# 
# 
# concatenated_sdata[table]=adata_vis
# 
# 
# concatenated_sdata[table] = mod.export_posterior(
#     concatenated_sdata[table], sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
# )
# 
# 
# 
# mod.plot_QC()
# plt.savefig('04-mod.plot_QC.png')
# 
# # print(adata_vis.obsm['q05_cell_abundance_w_sf'])
# # adata_vis.obsm['q05_cell_abundance_w_sf']
# # adata_vis.obsm['q95_cell_abundance_w_sf']
# 
# 
# concatenated_sdata[table].obs[concatenated_sdata[table].uns['mod']['factor_names']] = concatenated_sdata[table].obsm['q05_cell_abundance_w_sf']
# 
# print(concatenated_sdata)
# 
# 
# # select one slide
# from cell2location.utils import select_slide
# slide = select_slide(concatenated_sdata[table], 'Human_Lymph_Node')
# print("可用的 img_key 列表：", slide.uns['spatial'].keys())
# 
# 
# print(concatenated_sdata[table])
# print(concatenated_sdata[table].uns["mod"]["factor_names"])
# 
# 
# 
# 
# new_columns = concatenated_sdata[table].obs.columns.str.replace("+", "plus", regex=False)
# concatenated_sdata[table].obs.columns = new_columns
# print("\n修改后的 obs 列名：")
# print(concatenated_sdata[table].obs.columns.tolist())
# 
# new_col=[]
# for i in concatenated_sdata[table].uns["mod"]["factor_names"]:
#   new_col.append(i.replace("+", "plus"))
# 
# 
# 
# 
# 
# 
# existing_cell_cols = [col for col in new_col if col in concatenated_sdata[table].obs.columns]
# if len(existing_cell_cols) != len(new_col):
#     missing_cols = set(new_col) - set(existing_cell_cols)
#     print(f"警告：以下细胞类型列在 obs 中不存在，已自动过滤：{missing_cols}")
# 
# abundance_matrix = concatenated_sdata[table].obs[existing_cell_cols]
# 
# concatenated_sdata[table].obs["cellLoca_type"] = abundance_matrix.idxmax(axis=1)
# 
# # # （可选）同时计算每个 spot 的最大丰度值，便于后续分析
# # adata.obs["max_abundance_value"] = abundance_matrix.max(axis=1)
# 
# print("新增 cellLoca_type 列后的前 5 行数据：")
# print(concatenated_sdata[table].obs["cellLoca_type"].head())
# 
# 
# # concatenated_sdata.tables["table"] = adata
# 
# concatenated_sdata.pl.render_images(elements="Human_Lymph_Node_hires_image").pl.render_shapes(elements="Human_Lymph_Node", color='cellLoca_type').pl.show(coordinate_systems="Human_Lymph_Node_downscaled_hires")
# plt.savefig('sopoatial.png',dpi=300,bbox_inches='tight')
# 
# 
# 
# 
# 
# ############# each type 可视化
# 
# 
# concatenated_sdata[table].uns["spatial"]=adata.uns["spatial"]
# sc.pl.spatial(
#     concatenated_sdata[table],
#     color=existing_cell_cols,  # 细胞类型列名
#     show=False,img_key='hires',
#     save="_Cell2Loc_q05_cell_abundance_w_sf.png"  # 保存图表
# )
# 
# from cell2location.plt import plot_spatial
# clust_labels = new_col[:6]
# clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
# clust_col
# 
# slide = select_slide(adata_vis, 'Human_Lymph_Node')
# 
# with mpl.rc_context({'figure.figsize': (15, 15)}):
#     fig = plot_spatial(
#         adata=concatenated_sdata[table],
#         color=clust_col, labels=clust_labels,
#         show_img=True,
#         style='fast',
#         max_color_quantile=0.992,
#         circle_diameter=6,
#         colorbar_position='right'
#     )
# plt.savefig('06-sc.pl.spatial.png')





























# element = "Human_Lymph_Node"
# # concatenated_sdata = blobs_annotating_element(element)
# image = "Human_Lymph_Node_hires_image"
# affine = Affine([[1, 2, 3], [-2, 1, 6], [0, 0, 1]], input_axes=("x", "y"), output_axes=("x", "y"))
# set_transformation(concatenated_sdata[element], affine, "transformed")
# set_transformation(concatenated_sdata[image], affine, "transformed")
# adata_vis = to_legacy_anndata(concatenated_sdata, include_images=True, coordinate_system="transformed")
# 
# print(adata_vis.uns,"#################")



# concatenated_sdata.pl.render_images(elements="Human_Lymph_Node_hires_image").pl.render_shapes(elements="Human_Lymph_Node", color='B_naive').pl.show(coordinate_systems="Human_Lymph_Node_downscaled_hires")
# plt.savefig('sopoatial.png')

# with mpl.rc_context({'axes.facecolor':  'black','figure.figsize': [4.5, 5]}):
#     sc.pl.spatial(slide, cmap='magma',
#                   # show first 8 cell types
#                   color=sdata_st["table"].uns["mod"]["factor_names"],
#                   ncols=4, size=1.3,
#                   vmin=0, vmax='p99.2')
#                  
# 
# plt.savefig('05-sc.pl.spatial.png')
















































































# import os
# # this line should go before importing cell2location
# os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
# import torch
# from numpy.core.multiarray import _reconstruct
# import sys
# import scanpy as sc
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# import cell2location
# from matplotlib import rcParams
# rcParams['pdf.fonttype'] = 42
# 
# 
# 
# adata_vis = sc.read_visium("./data/V1_Human_Lymph_Node/", count_file='V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5', load_images=True)
# adata_vis
# adata_vis.var_names
# 
# # rename genes to ENSEMBL ID for correct matching between single cell and spatial data
# adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
# adata_vis.var['SYMBOL'] = adata_vis.var_names
# adata_vis.var.set_index('gene_ids', drop=True, inplace=True)
# adata_vis
# adata_vis.var_names
# 
# adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
# 
# # remove MT genes for spatial mapping (keeping their counts in the object)
# adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
# adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
# 
# 
# # Read scrna data
# adata_ref = sc.read('./data/V1_Human_Lymph_Node/sc.h5ad')
# 
# # 将数据的genes转为ENSEMBL ID，与空间数据对应同样的id
# adata_ref.var['SYMBOL'] = adata_ref.var.index
# # rename 'GeneID-2' as necessary for your data
# adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)
# adata_ref.var
# 
# 
# 
# from cell2location.utils.filtering import filter_genes
# selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# # filter the object
# adata_ref = adata_ref[:, selected].copy()
# 
# 
# 
# cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
#                         # 10X reaction / sample / batch
#                         batch_key='Sample',
#                         # cell type, covariate used for constructing signatures
#                         labels_key='Subset',
#                         # multiplicative technical effects (platform, 3' vs 5', donor effect)
#                         categorical_covariate_keys=['Method']
#                        )
# 
# # create the regression model
# from cell2location.models import RegressionModel
# mod = RegressionModel(adata_ref)
# 
# # view anndata_setup as a sanity check
# mod.view_anndata_setup()
# 
# mod.train(max_epochs=250)
# 
# 
# mod.plot_history(20)
# plt.savefig('01-mod.plot_history.png')
# 
# 
# adata_ref = mod.export_posterior(
#     adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
# )
# 
# # Save model
# mod.save('./reference_signatures/', overwrite=True)
# adata_file = "./reference_signatures/sc.h5ad"
# adata_ref.write(adata_file)
# adata_file
# 
# 
# adata_file = "./reference_signatures/sc.h5ad"
# adata_ref = sc.read_h5ad(adata_file)
# with torch.serialization.safe_globals([_reconstruct]):
#     mod = cell2location.models.RegressionModel.load(
#     './reference_signatures/', 
#     adata_ref,
#     load_kwargs={"weights_only": False}  # 传递参数给torch.load
# )
# mod.plot_QC()
# 
# 
# if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
#     inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
#                             for i in adata_ref.uns['mod']['factor_names']]].copy()
# else:
#     inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
#                             for i in adata_ref.uns['mod']['factor_names']]].copy()
#                             
# inf_aver.columns = adata_ref.uns['mod']['factor_names']
# inf_aver.iloc[0:12, 0:12]
# 
# 
# 
# 
# intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
# adata_vis = adata_vis[:, intersect].copy()
# inf_aver = inf_aver.loc[intersect, :].copy()
# 
# # prepare anndata for cell2location model
# cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
# 
# 
# 
# 
# # create and train the model
# mod = cell2location.models.Cell2location(
#     adata_vis, cell_state_df=inf_aver,
#     # the expected average cell abundance: tissue-dependent
#     # hyper-prior which can be estimated from paired histology:
#     N_cells_per_location=30,
#     # hyperparameter controlling normalisation of
#     # within-experiment variation in RNA detection:
#     detection_alpha=20
# )
# mod.view_anndata_setup()
# 
# 
# 
# 
# mod.train(max_epochs=30000,
#           # train using full data (batch_size=None)
#           batch_size=None,
#           # use all data points in training because
#           # we need to estimate cell abundance at all locations
#           train_size=1
#          )
# 
# # plot ELBO loss history during training, removing first 100 epochs from the plot
# mod.plot_history(1000)
# plt.legend(labels=['full data training']);
# plt.savefig('03-mod.plot_history.png')
# 
# 
# adata_vis = mod.export_posterior(
#     adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
# )
# 
# mod.save("./", overwrite=True)
# # 生成 model.pt
# 
# # Save anndata object with results
# adata_file = "./sp.h5ad"
# adata_vis.write(adata_file)
# adata_file
# 
# 
# 
# 
# adata_file = "./sp.h5ad"
# adata_vis = sc.read_h5ad(adata_file)
# mod = cell2location.models.Cell2location.load("./", adata_vis)
# 
# 
# mod.plot_QC()
# plt.savefig('04-mod.plot_QC.png')
# 
# 
# 
# 
# adata_vis.obsm
# adata_vis.obsm['q05_cell_abundance_w_sf']
# adata_vis.obsm['q95_cell_abundance_w_sf']
# 
# 
# adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# 
# # select one slide
# from cell2location.utils import select_slide
# slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')
# 
# # plot in spatial coordinates
# with mpl.rc_context({'axes.facecolor':  'black','figure.figsize': [4.5, 5]}):
#     sc.pl.spatial(slide, cmap='magma',
#                   # show first 8 cell types
#                   color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC','B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],
#                   ncols=4, size=1.3,img_key='hires',
#                   # limit color scale at 99.2% quantile of cell abundance
#                   vmin=0, vmax='p99.2'
#                  )
# plt.savefig('05-sc.pl.spatial.png') 
# 
# 
# 
# from cell2location.plt import plot_spatial
# 
# # select up to 6 clusters
# clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']
# clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
# clust_col
# 
# slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')
# 
# with mpl.rc_context({'figure.figsize': (15, 15)}):
#     fig = plot_spatial(
#         adata=slide,
#         # labels to show on a plot
#         color=clust_col, labels=clust_labels,
#         show_img=True,
#         # 'fast' (white background) or 'dark_background'
#         style='fast',
#         # limit color scale at 99.2% quantile of cell abundance
#         max_color_quantile=0.992,
#         # size of locations (adjust depending on figure size)
#         circle_diameter=6,
#         colorbar_position='right'
#     )
# plt.savefig('06-sc.pl.spatial.png')
