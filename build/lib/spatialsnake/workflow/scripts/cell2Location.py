import os
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from PIL import Image
from spatialdata.transformations import Identity, Scale
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
from matplotlib import use
from scipy.sparse import issparse
from scvi import REGISTRY_KEYS
from typing import Union, List, Optional, Iterable, Sequence, Dict
from matplotlib.axes import Axes
from anndata import AnnData
import matplotlib


def parse_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    value_str = str(value).strip().lower()
    if value_str in {"true", "1", "yes", "y", "t"}:
        return True
    if value_str in {"false", "0", "no", "n", "f", ""}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

def ensure_var_index(adata, candidates):
    for key in candidates:
        if key in adata.var.columns:
            adata.var.set_index(key, drop=True, inplace=True)
            return key
    return None

def ensure_symbol_column(adata):
    if "SYMBOL" not in adata.var.columns:
        adata.var["SYMBOL"] = adata.var_names


def load_spatial_input(path):
    if path.endswith(".h5ad"):
        adata = sc.read_h5ad(path)
        return adata, None, None
    sdata = spd.read_zarr(path)
    table_key = next(iter(sdata.tables.keys()))
    return sdata[table_key], sdata, table_key

def write_adata(adata, path):
    if path.endswith(".h5ad"):
        adata.write(path)
    else:
        adata.write_zarr(path)

parser = argparse.ArgumentParser()
parser.add_argument("--input_spatial", required=True)
parser.add_argument("--input_singlecell", required=True)
parser.add_argument("--output_dir_zarr", required=True)
parser.add_argument("--max_epochs_reference", type=int, default=250)
parser.add_argument("--remove_mt", type=parse_bool, default=True)
parser.add_argument("--sample_id", type=str, required=True)
parser.add_argument("--type", type=str, required=True)
parser.add_argument("--max_epochs_st", type=int, default=30000)
parser.add_argument("--N_cells_per_location", type=int, default=30)



parser.add_argument("--labels_key_reference", default="celltype", required=False)
parser.add_argument("--batch_key_reference", default="sample", required=False)
parser.add_argument("--cell_count_cutoff", default=15, required=False)
parser.add_argument("--cell_percentage_cutoff2", default=0.05, required=False)
parser.add_argument("--nonz_mean_cutoff", default=1.12, required=False)
parser.add_argument("--labels_key_st", default="celltype", required=False)
parser.add_argument("--batch_key_st", default="sample", required=False)
parser.add_argument("--layer_st", default=None, required=False)
parser.add_argument("--detection_alpha", type=float, default=20, required=False)
parser.add_argument("--device", required=False)
parser.add_argument("--save_models", type=parse_bool, default=True, required=False)
args = parser.parse_args()
if args.save_models:
    print("@@@@@@@@@@@@@@@@@@@555555")
if args.device == "cuda":
    os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
else:
    os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32'

def cell2loc_plot_history(model, fig_path, iter_start=0, iter_end=None, rolling_window=1):
# Adapted from: https://github.com/BayraktarLab/cell2location/blob/master/cell2location/models/base/_pyro_mixin.py#L407
    if not hasattr(model, "history_") or "elbo_train" not in model.history_:
        raise RuntimeError("Model history is unavailable, please run train() before plotting history.")

    train_history = model.history_["elbo_train"].copy()
    if iter_end is None:
        iter_end = len(train_history)

    train_history = train_history.iloc[iter_start:iter_end]
    if train_history.empty:
        raise ValueError("Training history is empty after applying iter_start/iter_end.")

    train_x = np.asarray(train_history.index, dtype=float)
    train_y = np.asarray(train_history.values).flatten().astype(float)

    if rolling_window and rolling_window > 1 and len(train_y) >= rolling_window:
        train_y = (
            pd.Series(train_y)
            .rolling(window=rolling_window, min_periods=1)
            .mean()
            .to_numpy()
        )

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(train_x, train_y, color="#1f77b4", linewidth=1.8, label="train")

    if "elbo_validation" in model.history_:
        validation_history = model.history_["elbo_validation"].copy().iloc[iter_start:iter_end]
        if not validation_history.empty:
            valid_x = np.asarray(validation_history.index, dtype=float)
            valid_y = np.asarray(validation_history.values).flatten().astype(float)
            if rolling_window and rolling_window > 1 and len(valid_y) >= rolling_window:
                valid_y = (
                    pd.Series(valid_y)
                    .rolling(window=rolling_window, min_periods=1)
                    .mean()
                    .to_numpy()
                )
            ax.plot(valid_x, valid_y, color="#ff7f0e", linewidth=1.8, label="validation")

    ax.set_xlim(train_x.min(), train_x.max() if train_x.max() > train_x.min() else train_x.min() + 1)
    ax.set_xlabel("Training epochs")
    ax.set_ylabel("-ELBO loss")
    ax.set_title("Training history")
    ax.grid(alpha=0.3, linestyle="--", linewidth=0.5)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

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

adata_vis, concatenated_sdata, table_key = load_spatial_input(args.input_spatial)

adata_vis.var['SYMBOL'] = adata_vis.var_names
ensure_var_index(adata_vis, ["gene_ids", "gene_id", "ensembl"])

def ensure_counts(adata, name="adata"):
    if issparse(adata.X):
        is_int = np.all(np.mod(adata.X.data, 1) == 0)
    else:
        is_int = np.all(np.mod(adata.X, 1) == 0)
    
    if not is_int:
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
        elif "raw_counts" in adata.layers:
            adata.X = adata.layers["raw_counts"].copy()
        elif adata.raw is not None:
             if adata.raw.shape == adata.shape:
                 adata.X = adata.raw.X.copy()
                 print(f"{name}.X  -> using {name}.raw.X")
    else:
        print(f"{name}.X  not integer, no need to ensure_counts")

ensure_counts(adata_vis, "adata_vis")

if "sample" not in adata_vis.obs:
    if "region" in adata_vis.obs:
        adata_vis.obs["sample"] = adata_vis.obs["region"].astype(str)
    elif "spatial" in adata_vis.uns and len(adata_vis.uns["spatial"]) > 0:
        adata_vis.obs["sample"] = list(adata_vis.uns["spatial"].keys())[0]
    else:
        adata_vis.obs["sample"] = str(args.sample_id)

ensure_symbol_column(adata_vis)
ensure_var_index(adata_vis, ["gene_ids", "gene_id", "ensembl"])

adata_ref = sc.read_h5ad(args.input_singlecell)
ensure_counts(adata_ref, "adata_ref")
ensure_symbol_column(adata_ref)
ensure_var_index(adata_ref, ["GeneID-2", "gene_ids", "gene_id", "ensembl"])
adata_ref.var['SYMBOL'] = adata_ref.var.index
print(adata_ref)
from cell2location.utils.filtering import filter_genes
if args.remove_mt:
    print("remove_mt")
    symbol_col = adata_vis.var.get("SYMBOL", adata_vis.var_names)
    mt_mask = symbol_col.str.startswith(("MT-", "mt-"))
    if mt_mask.any():
        adata_vis.var["MT_gene"] = mt_mask
        adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]
    shared_features = [f for f in adata_vis.var_names if f in adata_ref.var_names]
    adata_ref = adata_ref[:, shared_features].copy()
    adata_vis = adata_vis[:, shared_features].copy()
    selected = filter_genes(adata_ref,cell_count_cutoff=float(args.cell_count_cutoff),
    cell_percentage_cutoff2=float(args.cell_percentage_cutoff2),
    nonz_mean_cutoff=float(args.nonz_mean_cutoff))
    adata_ref = adata_ref[:, selected].copy()
    adata_vis = adata_vis[:, selected].copy()

print(adata_vis)
print(adata_ref)
print(adata_ref.X)
print(adata_ref.obs)


if args.labels_key_reference not in adata_ref.obs:
    raise ValueError(f"labels_key_reference not found in adata_ref.obs: {args.labels_key_reference}")
if args.batch_key_reference not in adata_ref.obs:
    raise ValueError(f"batch_key_reference not found in adata_ref.obs: {args.batch_key_reference}")
adata_ref.obs[args.labels_key_reference] = pd.Categorical(
    adata_ref.obs[args.labels_key_reference].astype(str).str.strip()
)
adata_ref.obs[args.batch_key_reference] = pd.Categorical(
    adata_ref.obs[args.batch_key_reference].astype(str).str.strip()
)
categorical_keys = ["Method"] if "Method" in adata_ref.obs else None
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    labels_key=args.labels_key_reference,
    batch_key=args.batch_key_reference,
    categorical_covariate_keys=categorical_keys
)

from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()






mod.train(max_epochs=args.max_epochs_reference)
cell2loc_plot_history(mod, output_dir + "/figure/ELBO_sc_model.png")

ref_batch_size = min(2500, adata_ref.n_obs)
adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": ref_batch_size
    }
)

if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]].copy()
else:
    inf_aver = adata_ref.var[[f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]].copy()
inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
inf_aver.to_csv(output_dir+"/Cell2Loc_inf_aver.csv")

print(adata_vis)
print(adata_ref)
print(adata_ref.X)
print(adata_ref.obs)



cell2loc_plot_QC_reference(mod, output_dir + "/figure/QC_reference_reconstruction_accuracy.png", output_dir + "/figure/QC_reference_expression signatures_vs_avg_expression.png")


if args.save_models:
  mod.save(output_dir +"/Reference_model", overwrite=True)



intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
batch_key_st = args.batch_key_st if args.batch_key_st in adata_vis.obs else None
if batch_key_st is not None:
    adata_vis.obs[batch_key_st] = pd.Categorical(adata_vis.obs[batch_key_st].astype(str).str.strip())
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=batch_key_st)



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


model_spatial.train(max_epochs=args.max_epochs_st, batch_size=None)
cell2loc_plot_history(model_spatial, output_dir + "/figure/ELBO_spatial_model.png")
spatial_batch_size = min(adata_vis.n_obs, model_spatial.adata.n_obs) if adata_vis.n_obs > 0 else 250
adata_vis = model_spatial.export_posterior(
    adata_vis, sample_kwargs={"num_samples": 1000, "batch_size": spatial_batch_size}
)




cell2loc_plot_QC_reconstr(model_spatial, output_dir + "/figure/QC_spatial_reconstruction_accuracy.png")

if args.save_models:
    model_spatial.save(output_dir +"/Spatial_model", overwrite=True)


if concatenated_sdata is not None:
    concatenated_sdata[table_key] = adata_vis
    concatenated_sdata.write(args.output_dir_zarr)
else:
    write_adata(adata_vis, args.output_dir_zarr)
