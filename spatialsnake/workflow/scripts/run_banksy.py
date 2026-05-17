import argparse
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as spd
import matplotlib.pyplot as plt
script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)
if workflow_dir not in sys.path:
    sys.path.insert(0, workflow_dir)
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.main import concatenate_all
from banksy_utils.umap_pca import pca_umap
from banksy.cluster_methods import run_Leiden_partition
from banksy.plot_banksy import plot_results
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, matthews_corrcoef
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.plot import plot_celltype_spatial_enrichment

parser = argparse.ArgumentParser(description="Run BANKSY clustering and save to spatialdata zarr")
parser.add_argument("--input_dir", type=str, required=True)
parser.add_argument("--output_zarr_path", type=str, required=True)
parser.add_argument("--k_geom", type=float, required=False)
parser.add_argument("--max_m", type=float, required=False)
parser.add_argument("--nbr_weight_decay", type=str, required=False)
parser.add_argument("--n_comps", type=str, required=False)
parser.add_argument("--lambda_list", type=str, required=False)
parser.add_argument("--resolution", type=str, required=False)
args = parser.parse_args()

def parse_list(value, cast=float, default=None):
    if value is None or value == "":
        return default
    if isinstance(value, list):
        return [cast(v) for v in value]
    text = str(value).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
    if "," in text:
        parts = [p for p in text.split(",") if p.strip() != ""]
        return [cast(p.strip()) for p in parts]
    return [cast(text)]

def ensure_spatial_coords(adata):
    if "spatial" not in adata.obsm:
        if "array_row" in adata.obs.columns and "array_col" in adata.obs.columns:
            adata.obsm["spatial"] = np.vstack(
                (adata.obs["array_row"].to_numpy(), adata.obs["array_col"].to_numpy())
            ).T
        else:
            raise ValueError("Missing spatial coordinates in adata")
    if "array_row" not in adata.obs.columns:
        adata.obs["array_row"] = adata.obsm["spatial"][:, 0]
    if "array_col" not in adata.obs.columns:
        adata.obs["array_col"] = adata.obsm["spatial"][:, 1]
    return ("array_row", "array_col", "spatial")

def load_spatial_input(path):
    if path.endswith(".h5ad"):
        return sc.read_h5ad(path), None, None
    sdata = spd.read_zarr(path)
    table_key = None
    adata = None
    for table in sdata.tables.keys():
        table_key = table
        adata = sdata[table]
        break
    if adata is None:
        raise ValueError("No tables found in spatialdata zarr")
    return adata, sdata, table_key

def write_output(adata, sdata, table_key, output_path):
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    if sdata is not None and table_key is not None and output_path.endswith(".zarr"):
        sdata[table_key] = adata
        sdata.write(output_path, overwrite=True)
    else:
        if output_path.endswith(".h5ad"):
            adata.write(output_path)
        else:
            adata.write_zarr(output_path)

k_geom = int(args.k_geom) if args.k_geom is not None else 15
max_m = int(args.max_m) if args.max_m is not None else 1
nbr_weight_decay = args.nbr_weight_decay if args.nbr_weight_decay else "scaled_gaussian"
pca_dims = parse_list(args.n_comps, cast=int, default=[20])
lambda_list = parse_list(args.lambda_list, cast=float, default=[0.8])
resolutions = parse_list(args.resolution, cast=float, default=[0.5])

adata, sdata, table_key = load_spatial_input(args.input_dir)
coord_keys = ensure_spatial_coords(adata)
output_dir = os.path.dirname(args.output_zarr_path) or "."
results_dir = os.path.join(output_dir, "banksy_results")
os.makedirs(results_dir, exist_ok=True)

banksy_dict = initialize_banksy(
    adata,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=False,
    plt_nbr_weights=False,
    plt_agf_angles=False,
    plt_theta=False
)

banksy_dict, _ = generate_banksy_matrix(
    adata,
    banksy_dict,
    lambda_list,
    max_m,
    plot_std=False,
    save_matrix=False,
    variance_balance=False,
    verbose=True
)

pca_umap(
    banksy_dict,
    pca_dims=pca_dims,
    add_umap=True,
    plt_remaining_var=False
)

results_df, max_num_labels = run_Leiden_partition(
    banksy_dict,
    resolutions,
    num_nn=50,
    num_iterations=-1,
    partition_seed=12345,
    match_labels=False
)

if results_df is None or results_df.shape[0] == 0:
    raise ValueError("BANKSY clustering did not return any results")

label_obj = results_df.iloc[0]["labels"]
labels = label_obj.dense if hasattr(label_obj, "dense") else np.asarray(label_obj)
if labels.shape[0] != adata.n_obs:
    raise ValueError("BANKSY labels length does not match number of observations")

adata.obs["spatial_cluster"] = pd.Categorical(labels)
results_df.to_csv(os.path.join(results_dir, "banksy_results.csv"))

if "celltype" in adata.obs.columns:
    plot_celltype_spatial_enrichment(
        adata,
        celltype_col="celltype",
        cluster_col="spatial_cluster",
        save_path=os.path.join(results_dir, "celltype_spatial_cluster_enrichment.png"),
        title="Cell type enrichment across BANKSY spatial clusters"
    )

decay_key = nbr_weight_decay if nbr_weight_decay in banksy_dict else next(iter(banksy_dict.keys()))
weights_by_m = banksy_dict[decay_key]["weights"]
weights_graph = weights_by_m[max_m] if max_m in weights_by_m else weights_by_m[0]
plot_results(
    results_df,
    weights_graph,
    "tab20",
    match_labels=False,
    coord_keys=coord_keys,
    max_num_labels=max_num_labels,
    save_path=os.path.join(results_dir, "BANKSY-Results"),
    save_fig=False,
    save_fullfig=True,
    dataset_name="BANKSY",
    save_labels=True
)

if "celltype" in adata.obs.columns:
    cmap = plt.get_cmap("tab20")
    sc.pl.scatter(
        adata,
        x="array_row",
        y="array_col",
        color="celltype",
        color_map=cmap,
        title="Tissue",
        size=5,
        show=False
    )
    plt.savefig(os.path.join(results_dir, "scatter.png"), dpi=300, bbox_inches="tight")
    plt.close()

    nonspatial_dict = {"nonspatial": {0.0: {"adata": concatenate_all([adata.X], 0, adata=adata)}}}
    pca_umap(nonspatial_dict, pca_dims=pca_dims, add_umap=True)
    nonspatial_df, nonspatial_max_num_labels = run_Leiden_partition(
        nonspatial_dict,
        resolutions,
        num_nn=50,
        num_iterations=-1,
        partition_seed=12345,
        match_labels=False
    )
    plot_results(
        nonspatial_df,
        weights_graph,
        "tab20",
        match_labels=False,
        coord_keys=coord_keys,
        max_num_labels=nonspatial_max_num_labels,
        save_path=os.path.join(results_dir, "BANKSY-Results-Nonspatial"),
        save_fig=False,
        save_fullfig=True,
        dataset_name="BANKSY-Nonspatial",
        save_labels=True
    )

    banksy_spatial_clusters = results_df.iloc[0]["labels"]
    nonspatial_clusters = nonspatial_df.iloc[0]["labels"]

    def calculate_metrics(cluster_labels, annotated_labels):
        ari_score = adjusted_rand_score(cluster_labels, annotated_labels)
        ami_score = adjusted_mutual_info_score(cluster_labels, annotated_labels)
        if isinstance(annotated_labels.dtype, pd.CategoricalDtype):
            annotated_labels = annotated_labels.cat.codes
        mcc_score = matthews_corrcoef(cluster_labels, annotated_labels)
        return ari_score, ami_score, mcc_score

    nonspatial_ari, nonspatial_ami, nonspatial_mcc = calculate_metrics(
        nonspatial_clusters.dense, adata.obs["celltype"]
    )
    banksy_ari, banksy_ami, banksy_mcc = calculate_metrics(
        banksy_spatial_clusters.dense, adata.obs["celltype"]
    )

    def bar_plot(metrics, methods):
        fig, ax = plt.subplots(figsize=(8, 8), layout="constrained")
        x = np.arange(len(methods))
        width = 0.25
        multiplier = 0
        for method, metric in metrics.items():
            offset = width * multiplier
            rects = ax.bar(x + offset, metric, width, label=method)
            ax.bar_label(rects, padding=3)
            multiplier += 1
        ax.set_ylabel("Metrices", fontsize=18)
        ax.set_title("Similarity between BANKSY labels and annotated communities", fontsize=20)
        ax.set_xticks(x + width, methods, fontsize=16)
        ax.legend(loc="upper left", fontsize=18)
        plt.savefig(os.path.join(results_dir, "bar.png"), dpi=300, bbox_inches="tight")
        plt.close()

    metrics = {
        "Adjusted Rand Index": (nonspatial_ari, banksy_ari),
        "Adjusted Mutual Information": (nonspatial_ami, banksy_ami),
        "Matthew Correlation Coefficient": (nonspatial_mcc, banksy_mcc),
    }
    bar_plot(metrics, ("Non-spatial labels", "BANKSY labels"))

sample_id = os.path.splitext(os.path.basename(args.output_zarr_path))[0]
data_type = "xenium" if "xenium" in f"{args.input_dir}{args.output_zarr_path}".lower() else "visium"
export_cluster_csv(
    sdata if sdata is not None else adata,
    data_type,
    output_dir,
    cell_id_col="cell_id",
    info_col="spatial_cluster",
    sample_col="sample",
    sample_id=sample_id
)

write_output(adata, sdata, table_key, args.output_zarr_path)
