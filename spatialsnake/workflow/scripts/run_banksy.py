import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as spd
from scipy import sparse

script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)
if workflow_dir not in sys.path:
    sys.path.insert(0, workflow_dir)

from banksy.cluster_methods import run_Leiden_partition
from banksy.embed_banksy import generate_banksy_matrix
from banksy.initialize_banksy import initialize_banksy
from banksy.main import concatenate_all
from banksy.plot_banksy import plot_results
from banksy_utils.umap_pca import pca_umap
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, matthews_corrcoef

from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.logging_utils import log_step, setup_logger
from spatialsnake.workflow.function.plot import plot_celltype_spatial_enrichment

logger = setup_logger("banksy")


parser = argparse.ArgumentParser(description="Run lightweight BANKSY clustering and save SpatialData/AnnData output")
parser.add_argument("--input_dir", type=str, required=True)
parser.add_argument("--output_zarr_path", type=str, required=True)
parser.add_argument("--k_geom", type=float, required=False)
parser.add_argument("--max_m", type=float, required=False)
parser.add_argument("--nbr_weight_decay", type=str, required=False)
parser.add_argument("--n_comps", type=str, required=False)
parser.add_argument("--resolution", type=str, required=False)
parser.add_argument("--lambda_list", type=str, required=False)
parser.add_argument("--banksy_n_comps", type=str, required=False)
parser.add_argument("--banksy_resolution", type=str, required=False)
parser.add_argument("--banksy_num_nn", type=int, required=False)
parser.add_argument("--banksy_max_features", type=int, required=False)
parser.add_argument("--banksy_feature_col", type=str, required=False)
parser.add_argument("--banksy_add_umap", type=str, required=False)
parser.add_argument("--banksy_plot_full", type=str, required=False)
parser.add_argument("--banksy_run_nonspatial", type=str, required=False)
parser.add_argument("--banksy_plot_celltype_enrichment", type=str, required=False)
parser.add_argument("--banksy_plot_max_points", type=int, required=False)
parser.add_argument("--banksy_sample_col", type=str, required=False)
parser.add_argument("--banksy_selected_lambda", type=str, required=False)
parser.add_argument("--banksy_selected_resolution", type=str, required=False)
parser.add_argument("--banksy_seed", type=int, required=False)
args = parser.parse_args()


def parse_bool(value, default=False):
    if value is None or value == "":
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y", "t"}:
        return True
    if text in {"false", "0", "no", "n", "f"}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def parse_list(value, cast=float, default=None):
    if value is None or value == "":
        return default
    if isinstance(value, (list, tuple)):
        return [cast(v) for v in value]
    text = str(value).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
    text = text.replace(";", ",")
    if "," in text:
        parts = [p for p in text.split(",") if p.strip() != ""]
        return [cast(p.strip()) for p in parts]
    return [cast(text)]


def parse_optional_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"none", "null", "na"}:
        return None
    return float(text)


def ensure_spatial_coords(adata):
    if "spatial" not in adata.obsm:
        if "array_row" in adata.obs.columns and "array_col" in adata.obs.columns:
            adata.obsm["spatial"] = np.vstack(
                (adata.obs["array_row"].to_numpy(), adata.obs["array_col"].to_numpy())
            ).T
        else:
            raise ValueError("Missing spatial coordinates in adata.obsm['spatial'] or obs array_row/array_col")
    if "array_row" not in adata.obs.columns:
        adata.obs["array_row"] = np.asarray(adata.obsm["spatial"])[:, 0]
    if "array_col" not in adata.obs.columns:
        adata.obs["array_col"] = np.asarray(adata.obsm["spatial"])[:, 1]
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
    elif output_path.endswith(".h5ad"):
        adata.write(output_path)
    else:
        adata.write_zarr(output_path)


def matrix_column_variance(x):
    if sparse.issparse(x):
        mean = np.asarray(x.mean(axis=0)).ravel()
        squared = x.copy()
        squared.data **= 2
        mean_sq = np.asarray(squared.mean(axis=0)).ravel()
        return mean_sq - mean**2
    if hasattr(x, "compute"):
        x = x.compute()
    arr = np.asarray(x)
    return np.var(arr, axis=0)


def select_feature_names(adata, max_features, feature_col):
    max_features = int(max_features)
    if max_features <= 0 or adata.n_vars <= max_features:
        return adata.var_names.to_numpy(), "all_features"

    if feature_col and feature_col in adata.var.columns:
        mask = adata.var[feature_col].astype(bool).to_numpy()
        selected = adata.var_names[mask].to_numpy()
        if selected.size > 0:
            if selected.size > max_features:
                rank_col = f"{feature_col}_rank"
                if rank_col in adata.var.columns:
                    ranks = pd.to_numeric(adata.var.loc[selected, rank_col], errors="coerce")
                    selected = ranks.sort_values(kind="mergesort").index.to_numpy()[:max_features]
                else:
                    selected = selected[:max_features]
            return selected, feature_col

    variances = matrix_column_variance(adata.X)
    if variances.size != adata.n_vars:
        raise ValueError("Could not compute per-gene variance for BANKSY feature selection")
    top_idx = np.argsort(np.nan_to_num(variances, nan=-np.inf))[-max_features:]
    top_idx = top_idx[np.argsort(np.nan_to_num(variances[top_idx], nan=-np.inf))[::-1]]
    return adata.var_names[top_idx].to_numpy(), "variance_top"


def label_to_array(label_obj):
    labels = label_obj.dense if hasattr(label_obj, "dense") else np.asarray(label_obj)
    return np.asarray(labels).reshape(-1)


def tidy_results(results_df, sample_name, selected_key=None):
    drop_cols = [col for col in ["labels", "adata"] if col in results_df.columns]
    summary = results_df.drop(columns=drop_cols).copy()
    summary.insert(0, "result_key", summary.index.astype(str))
    summary.insert(0, "sample", sample_name)
    summary["selected"] = summary["result_key"].astype(str) == str(selected_key)
    return summary.reset_index(drop=True)


def select_result_row(results_df, selected_lambda, selected_resolution, sample_name):
    mask = pd.Series(True, index=results_df.index)
    if selected_lambda is not None:
        if "lambda_param" not in results_df.columns:
            raise ValueError("BANKSY result table does not contain lambda_param")
        mask &= np.isclose(results_df["lambda_param"].astype(float), float(selected_lambda))
    if selected_resolution is not None:
        if "resolution" not in results_df.columns:
            raise ValueError("BANKSY result table does not contain resolution")
        mask &= np.isclose(results_df["resolution"].astype(float), float(selected_resolution))
    matches = results_df.loc[mask]
    if matches.empty:
        raise ValueError(
            f"No BANKSY result matched selected lambda={selected_lambda} "
            f"and resolution={selected_resolution} for sample {sample_name!r}"
        )
    selected_key = matches.index[0]
    return selected_key, matches.iloc[0]


def estimate_dense_banksy_gb(n_obs, n_features, max_m):
    banksy_blocks = 1 + max_m + 1
    return n_obs * n_features * banksy_blocks * 8 / (1024**3)


def run_banksy_single(adata_input, sample_name, feature_names, output_prefix):
    if adata_input.n_obs <= 2:
        raise ValueError(f"BANKSY sample {sample_name!r} has too few observations: {adata_input.n_obs}")
    local_k_geom = min(k_geom, adata_input.n_obs - 1)
    local_num_nn = min(banksy_num_nn, adata_input.n_obs - 1)
    adata_banksy = adata_input[:, feature_names].copy()
    ensure_spatial_coords(adata_banksy)

    est_gb = estimate_dense_banksy_gb(adata_banksy.n_obs, adata_banksy.n_vars, max_m)
    logger.info(
        f"BANKSY sample={sample_name}: n_obs={adata_banksy.n_obs}, "
        f"n_features={adata_banksy.n_vars}, k_geom={local_k_geom}, "
        f"estimated dense BANKSY matrix={est_gb:.2f} GB"
    )

    banksy_dict = initialize_banksy(
        adata_banksy,
        coord_keys,
        local_k_geom,
        nbr_weight_decay=nbr_weight_decay,
        max_m=max_m,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False,
    )

    banksy_dict, _ = generate_banksy_matrix(
        adata_banksy,
        banksy_dict,
        lambda_list,
        max_m,
        plot_std=False,
        save_matrix=False,
        variance_balance=False,
        verbose=False,
    )

    pca_umap(
        banksy_dict,
        pca_dims=pca_dims,
        add_umap=banksy_add_umap,
        plt_remaining_var=False,
    )

    results_df, max_num_labels = run_Leiden_partition(
        banksy_dict,
        resolutions,
        num_nn=local_num_nn,
        num_iterations=-1,
        partition_seed=banksy_seed,
        match_labels=False,
    )
    if results_df is None or results_df.shape[0] == 0:
        raise ValueError(f"BANKSY clustering did not return any results for sample {sample_name!r}")

    selected_key, selected_row = select_result_row(
        results_df,
        banksy_selected_lambda,
        banksy_selected_resolution,
        sample_name,
    )
    labels = label_to_array(selected_row["labels"])
    if labels.shape[0] != adata_input.n_obs:
        raise ValueError(
            f"BANKSY labels length does not match observations for sample {sample_name!r}: "
            f"{labels.shape[0]} vs {adata_input.n_obs}"
        )

    if banksy_plot_full:
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
            save_path=output_prefix,
            save_fig=False,
            save_fullfig=True,
            dataset_name="BANKSY",
            save_labels=True,
        )
        if banksy_run_nonspatial:
            run_nonspatial_baseline(adata_banksy, weights_graph, output_prefix, sample_name, labels)

    return labels, tidy_results(results_df, sample_name, selected_key)


def run_nonspatial_baseline(adata_banksy, weights_graph, output_prefix, sample_name, banksy_labels):
    if "celltype" not in adata_banksy.obs.columns:
        logger.info("Skipping non-spatial baseline metrics because obs['celltype'] is absent")
        return
    nonspatial_dict = {"nonspatial": {0.0: {"adata": concatenate_all([adata_banksy.X], 0, adata=adata_banksy)}}}
    pca_umap(nonspatial_dict, pca_dims=pca_dims, add_umap=banksy_add_umap, plt_remaining_var=False)
    nonspatial_df, nonspatial_max_num_labels = run_Leiden_partition(
        nonspatial_dict,
        resolutions,
        num_nn=min(banksy_num_nn, adata_banksy.n_obs - 1),
        num_iterations=-1,
        partition_seed=banksy_seed,
        match_labels=False,
    )
    plot_results(
        nonspatial_df,
        weights_graph,
        "tab20",
        match_labels=False,
        coord_keys=coord_keys,
        max_num_labels=nonspatial_max_num_labels,
        save_path=f"{output_prefix}-Nonspatial",
        save_fig=False,
        save_fullfig=True,
        dataset_name="BANKSY-Nonspatial",
        save_labels=True,
    )
    nonspatial_labels = label_to_array(nonspatial_df.iloc[0]["labels"])
    metrics = []
    for method, labels in [("Non-spatial labels", nonspatial_labels), ("BANKSY labels", banksy_labels)]:
        annotated = adata_banksy.obs["celltype"]
        annotated_codes = annotated.cat.codes if isinstance(annotated.dtype, pd.CategoricalDtype) else annotated
        metrics.append(
            {
                "sample": sample_name,
                "method": method,
                "adjusted_rand_index": adjusted_rand_score(labels, annotated),
                "adjusted_mutual_information": adjusted_mutual_info_score(labels, annotated),
                "matthews_correlation_coefficient": matthews_corrcoef(labels, annotated_codes),
            }
        )
    metrics_path = os.path.join(results_dir, "banksy_nonspatial_metrics.csv")
    metrics_df = pd.DataFrame(metrics)
    if os.path.exists(metrics_path):
        metrics_df.to_csv(metrics_path, mode="a", header=False, index=False)
    else:
        metrics_df.to_csv(metrics_path, index=False)


def get_sample_groups(adata, sample_col):
    if sample_col and sample_col in adata.obs.columns:
        values = adata.obs[sample_col].astype(str)
        unique = [value for value in pd.unique(values) if value not in {"", "nan", "None"}]
        if len(unique) > 1:
            return [(value, values == value) for value in unique], True
    return [("all", pd.Series(True, index=adata.obs_names))], False


def plot_spatial_clusters(adata, cluster_col, path, max_points):
    df = pd.DataFrame(
        {
            "x": np.asarray(adata.obs["array_row"], dtype=float),
            "y": np.asarray(adata.obs["array_col"], dtype=float),
            "cluster": adata.obs[cluster_col].astype(str).to_numpy(),
        },
        index=adata.obs_names,
    )
    if max_points and max_points > 0 and df.shape[0] > max_points:
        df_plot = df.sample(n=max_points, random_state=banksy_seed)
        title_suffix = f" (downsampled to {max_points:,} points)"
    else:
        df_plot = df
        title_suffix = ""
    categories = sorted(pd.unique(df_plot["cluster"]), key=lambda x: str(x))
    cmap = plt.get_cmap("tab20", max(1, len(categories)))
    color_map = {cat: cmap(i % cmap.N) for i, cat in enumerate(categories)}
    fig_width = 8
    fig_height = 8
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    point_size = max(0.2, min(5.0, 50000 / max(1, df_plot.shape[0])))
    for cat in categories:
        part = df_plot[df_plot["cluster"] == cat]
        ax.scatter(part["x"], part["y"], s=point_size, color=color_map[cat], label=cat, linewidths=0)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_title(f"BANKSY spatial clusters{title_suffix}")
    ax.set_xlabel("array_row")
    ax.set_ylabel("array_col")
    if len(categories) <= 30:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=4, frameon=False, fontsize=7)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_parameters(path, feature_source, feature_count, is_multisample):
    rows = [
        ("k_geom", k_geom),
        ("max_m", max_m),
        ("nbr_weight_decay", nbr_weight_decay),
        ("lambda_list", ",".join(map(str, lambda_list))),
        ("banksy_n_comps", ",".join(map(str, pca_dims))),
        ("banksy_resolution", ",".join(map(str, resolutions))),
        ("banksy_num_nn", banksy_num_nn),
        ("banksy_max_features", banksy_max_features),
        ("banksy_feature_col", banksy_feature_col),
        ("selected_feature_source", feature_source),
        ("selected_feature_count", feature_count),
        ("banksy_add_umap", banksy_add_umap),
        ("banksy_plot_full", banksy_plot_full),
        ("banksy_run_nonspatial", banksy_run_nonspatial),
        ("banksy_plot_celltype_enrichment", banksy_plot_celltype_enrichment),
        ("banksy_plot_max_points", banksy_plot_max_points),
        ("banksy_sample_col", banksy_sample_col),
        ("multi_sample_mode", is_multisample),
        ("banksy_selected_lambda", banksy_selected_lambda if banksy_selected_lambda is not None else ""),
        ("banksy_selected_resolution", banksy_selected_resolution if banksy_selected_resolution is not None else ""),
        ("banksy_seed", banksy_seed),
    ]
    pd.DataFrame(rows, columns=["parameter", "value"]).to_csv(path, index=False)


k_geom = int(args.k_geom) if args.k_geom is not None else 15
max_m = int(args.max_m) if args.max_m is not None else 1
nbr_weight_decay = args.nbr_weight_decay if args.nbr_weight_decay else "scaled_gaussian"
pca_dims = parse_list(args.banksy_n_comps or args.n_comps, cast=int, default=[20])
lambda_list = parse_list(args.lambda_list, cast=float, default=[0.8])
resolutions = parse_list(args.banksy_resolution or args.resolution, cast=float, default=[0.5])
banksy_num_nn = int(args.banksy_num_nn) if args.banksy_num_nn is not None else 50
banksy_max_features = int(args.banksy_max_features) if args.banksy_max_features is not None else 2000
banksy_feature_col = args.banksy_feature_col if args.banksy_feature_col else "highly_variable"
banksy_add_umap = parse_bool(args.banksy_add_umap, default=False)
banksy_plot_full = parse_bool(args.banksy_plot_full, default=False)
banksy_run_nonspatial = parse_bool(args.banksy_run_nonspatial, default=False)
banksy_plot_celltype_enrichment = parse_bool(args.banksy_plot_celltype_enrichment, default=True)
banksy_plot_max_points = int(args.banksy_plot_max_points) if args.banksy_plot_max_points is not None else 200000
banksy_sample_col = args.banksy_sample_col if args.banksy_sample_col else "region"
banksy_selected_lambda = parse_optional_float(args.banksy_selected_lambda)
banksy_selected_resolution = parse_optional_float(args.banksy_selected_resolution)
banksy_seed = int(args.banksy_seed) if args.banksy_seed is not None else 12345

np.random.seed(banksy_seed)

log_step(logger, 1, 5, "loading spatial input")
adata, sdata, table_key = load_spatial_input(args.input_dir)
coord_keys = ensure_spatial_coords(adata)
output_dir = os.path.dirname(args.output_zarr_path) or "."
results_dir = os.path.join(output_dir, "banksy_results")
os.makedirs(results_dir, exist_ok=True)
logger.info(f"Loaded {adata.n_obs} observations and {adata.n_vars} genes")
if "spatial_cluster" in adata.obs.columns:
    logger.info("obs['spatial_cluster'] already exists and will be overwritten by BANKSY labels")

feature_names, feature_source = select_feature_names(adata, banksy_max_features, banksy_feature_col)
logger.info(f"Selected {len(feature_names)} genes for BANKSY using source={feature_source}")

log_step(logger, 2, 5, "running BANKSY spatial graph, matrix, PCA and Leiden clustering")
sample_groups, is_multisample = get_sample_groups(adata, banksy_sample_col)
all_labels = pd.Series(index=adata.obs_names, dtype="object")
summary_frames = []

for sample_name, mask in sample_groups:
    obs_names = adata.obs_names[mask.to_numpy() if hasattr(mask, "to_numpy") else mask]
    output_prefix = os.path.join(results_dir, f"BANKSY-Results-{sample_name}" if is_multisample else "BANKSY-Results")
    labels, summary = run_banksy_single(adata[obs_names].copy(), sample_name, feature_names, output_prefix)
    labels = labels.astype(str)
    if is_multisample:
        labels = np.array([f"{sample_name}_{label}" for label in labels], dtype=object)
    all_labels.loc[obs_names] = labels
    summary_frames.append(summary)

if all_labels.isna().any():
    missing = int(all_labels.isna().sum())
    raise ValueError(f"BANKSY did not assign labels for {missing} observations")

adata.obs["spatial_cluster"] = pd.Categorical(all_labels.astype(str))
results_summary = pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame()
results_summary.to_csv(os.path.join(results_dir, "banksy_results.csv"), index=False)

log_step(logger, 3, 5, "saving BANKSY summary tables")
cluster_assignments = pd.DataFrame(
    {
        "obs_id": adata.obs_names,
        "spatial_cluster": adata.obs["spatial_cluster"].astype(str).to_numpy(),
    }
)
if banksy_sample_col in adata.obs.columns:
    cluster_assignments.insert(1, banksy_sample_col, adata.obs[banksy_sample_col].astype(str).to_numpy())
cluster_assignments.to_csv(os.path.join(results_dir, "banksy_cluster_assignments.csv"), index=False)
write_parameters(
    os.path.join(results_dir, "banksy_parameters.csv"),
    feature_source=feature_source,
    feature_count=len(feature_names),
    is_multisample=is_multisample,
)

log_step(logger, 4, 5, "saving lightweight BANKSY visualizations")
plot_spatial_clusters(
    adata,
    "spatial_cluster",
    os.path.join(results_dir, "banksy_spatial_cluster.png"),
    banksy_plot_max_points,
)
if banksy_plot_celltype_enrichment and "celltype" in adata.obs.columns:
    plot_celltype_spatial_enrichment(
        adata,
        celltype_col="celltype",
        cluster_col="spatial_cluster",
        save_path=os.path.join(results_dir, "celltype_spatial_cluster_enrichment.png"),
        title="Cell type enrichment across BANKSY spatial clusters",
    )

sample_id = os.path.splitext(os.path.basename(args.output_zarr_path))[0]
data_type = "xenium" if "xenium" in f"{args.input_dir}{args.output_zarr_path}".lower() else "visium"
export_cluster_csv(
    sdata if sdata is not None else adata,
    data_type,
    output_dir,
    cell_id_col="cell_id",
    info_col="spatial_cluster",
    sample_col="sample",
    sample_id=sample_id,
)

log_step(logger, 5, 5, f"saving BANKSY output to {args.output_zarr_path}")
write_output(adata, sdata, table_key, args.output_zarr_path)
logger.info("BANKSY module completed")
