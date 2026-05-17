import os
import sys
import inspect
import argparse
import spatialdata as spd
import scanpy as sc
import pandas as pd
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir in sys.path:
    sys.path.remove(script_dir)
import liana as li
from matplotlib import pyplot as plt

plt.rcParams["figure.dpi"] = 50

parser = argparse.ArgumentParser(description="Process spatial data and convert to zarr format")
parser.add_argument("--input_dir", type=str, required=True)
parser.add_argument("--sample_id", type=str, required=False)
parser.add_argument("--output_zarr_path", type=str, required=True)
parser.add_argument("--type", type=str, required=True)
parser.add_argument("--method", type=str, required=False, default="cellphonedb")
parser.add_argument("--resource_name", type=str, required=False, default="consensus")
parser.add_argument("--celltype", type=str, required=False, default="celltype")
parser.add_argument("--expr_prop", type=float, required=False, default=0.1)
parser.add_argument("--min_cells", type=int, required=False, default=5)
parser.add_argument("--use_raw", type=str, required=False, default="true")
parser.add_argument("--sp_cellchat", type=str, required=False)
args = parser.parse_args()
output_dir = os.path.dirname(args.output_zarr_path)
if output_dir == "":
    output_dir = "."
os.makedirs(output_dir, exist_ok=True)

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

def get_method_config(method_name, res_df):
    configs = {
        "cellphonedb": {"magnitude": "lr_means", "specificity": "cellphone_pvals"},
        "connectome": {"magnitude": "expr_prod", "specificity": "scaled_weight"},
        "logfc": {"magnitude": "lr_logfc", "specificity": None},
        "natmi": {"magnitude": "expr_prod", "specificity": "spec_weight"},
        "singlecellsignalr": {"magnitude": "lrscore", "specificity": None},
        "rank_aggregate": {"magnitude": "magnitude_rank", "specificity": "specificity_rank"},
        "geometric_mean": {"magnitude": "lr_gmeans", "specificity": "gmean_pvals"},
        "scseqcomm": {"magnitude": "inter_score", "specificity": None},
        "cellchat": {"magnitude": "lr_probs", "specificity": "cellchat_pvals"},
    }
    config = configs.get(method_name, {})
    magnitude = config.get("magnitude")
    specificity = config.get("specificity")
    if res_df is not None:
        if magnitude is None or magnitude not in res_df.columns:
            for candidate in ["lr_means", "expr_prod", "lr_logfc", "lr_gmeans", "lr_probs", "inter_score", "magnitude_rank", "lrscore"]:
                if candidate in res_df.columns:
                    magnitude = candidate
                    break
        if specificity is None or specificity not in res_df.columns:
            for candidate in ["cellphone_pvals", "scaled_weight", "spec_weight", "gmean_pvals", "cellchat_pvals", "specificity_rank"]:
                if candidate in res_df.columns:
                    specificity = candidate
                    break
    return magnitude, specificity

def run_method(adata, method_name, uns_key):
    method = getattr(li.mt, method_name, None)
    if method is None:
        available = [m for m in dir(li.mt) if not m.startswith("_")]
        raise ValueError(f"Unsupported method: {method_name}. Available methods: {available}")
    use_raw = parse_bool(args.use_raw)
    if use_raw and adata.raw is None:
        print("adata.raw is not initialized, switching use_raw to False")
        use_raw = False
    args_map = {
        "adata": adata,
        "groupby": args.celltype,
        "resource_name": args.resource_name,
        "expr_prop": args.expr_prop,
        "min_cells": args.min_cells,
        "use_raw": use_raw,
        "key_added": uns_key,
        "verbose": True,
    }
    sig = inspect.signature(method)
    call_args = {k: v for k, v in args_map.items() if k in sig.parameters}
    method(**call_args)

def plot_results(adata, method_name, uns_key, output_dir):
    res_df = adata.uns.get(uns_key)
    magnitude, specificity = get_method_config(method_name, res_df)
    if magnitude is None:
        return
    adata.obs[args.celltype] = pd.Categorical(adata.obs[args.celltype].astype(str).str.strip())
    celltypes = adata.obs[args.celltype].cat.categories.astype(str).tolist()
    if len(celltypes) == 0:
        return
    plot_uns_key = uns_key
    if isinstance(res_df, pd.DataFrame) and "source" in res_df.columns and "target" in res_df.columns:
        valid_cells = set(celltypes)
        res_df_plot = res_df.copy()
        res_df_plot["source"] = res_df_plot["source"].astype(str).str.strip()
        res_df_plot["target"] = res_df_plot["target"].astype(str).str.strip()
        res_df_plot = res_df_plot[res_df_plot["source"].isin(valid_cells) & res_df_plot["target"].isin(valid_cells)]
        if res_df_plot.shape[0] == 0:
            print("Skip LIANA plotting: no source/target pairs matched current celltype categories")
            return
        plot_uns_key = f"{uns_key}__plot"
        adata.uns[plot_uns_key] = res_df_plot
        res_df = res_df_plot
    source_labels = celltypes[:2] if len(celltypes) > 1 else celltypes
    target_labels = celltypes[-2:] if len(celltypes) > 1 else celltypes
    if isinstance(res_df, pd.DataFrame):
        if "source" in res_df.columns:
            available_sources = set(res_df["source"].astype(str).unique().tolist())
            source_labels = [s for s in source_labels if s in available_sources]
        if "target" in res_df.columns:
            available_targets = set(res_df["target"].astype(str).unique().tolist())
            target_labels = [t for t in target_labels if t in available_targets]
        if len(source_labels) == 0 or len(target_labels) == 0:
            return
    filter_fun = None
    inverse_size = False
    orderby_ascending = False
    inverse_score = False
    if specificity and "pval" in specificity:
        def filter_fun(x):
            return x[specificity] <= 0.05
        inverse_size = True
        orderby_ascending = True
        if method_name == "cellphonedb":
            inverse_score = True
    size_key = specificity if specificity else magnitude
    dotplot_fig = li.pl.dotplot(
        adata=adata,
        colour=magnitude,
        size=size_key,
        inverse_size=inverse_size,
        source_labels=source_labels,
        target_labels=target_labels,
        figure_size=(12, 10),
        filter_fun=filter_fun,
        uns_key=plot_uns_key
    )
    dotplot_fig.save(os.path.join(output_dir, "dotplot.png"), dpi=300)
    my_plot = li.pl.tileplot(
        adata=adata,
        fill=magnitude,
        label=size_key,
        label_fun=lambda x: f"{x:.2f}",
        top_n=10,
        orderby=size_key,
        orderby_ascending=orderby_ascending,
        source_labels=source_labels,
        target_labels=target_labels,
        uns_key=plot_uns_key,
        source_title="Ligand",
        target_title="Receptor",
        figure_size=(10, 8)
    )
    my_plot.save(os.path.join(output_dir, "tileplot.png"), dpi=300)
    try:
        li.pl.circle_plot(
            adata,
            groupby=args.celltype,
            score_key=magnitude,
            inverse_score=inverse_score,
            source_labels=source_labels,
            filter_fun=filter_fun,
            pivot_mode="counts",
            figure_size=(10, 10),
            uns_key=plot_uns_key
        )
    except KeyError:
        try:
            li.pl.circle_plot(
                adata,
                groupby=args.celltype,
                score_key=magnitude,
                inverse_score=inverse_score,
                source_labels=celltypes,
                filter_fun=None,
                pivot_mode="counts",
                figure_size=(10, 10),
                uns_key=plot_uns_key
            )
        except KeyError:
            print("Skip circle plot due to unmatched source labels after fallback")
            return
    plt.savefig(os.path.join(output_dir, "circle.png"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()

method_aliases = {
    "log2fc": "logfc",
    "cellphone_db": "cellphonedb"
}
method_name = args.method.lower()
if method_name in method_aliases:
    method_name = method_aliases[method_name]
uns_key = f"{method_name}_res"
input_ext = os.path.splitext(args.input_dir)[1].lower()
input_is_h5ad = input_ext == ".h5ad"
table_key = None
concatenated_sdata = None

if input_is_h5ad:
    adata = sc.read_h5ad(args.input_dir)
else:
    concatenated_sdata = spd.read_zarr(args.input_dir)
    table_key = next(iter(concatenated_sdata.tables.keys()))
    adata = concatenated_sdata[table_key]

run_method(adata, method_name, uns_key)
plot_results(adata, method_name, uns_key, output_dir)

if input_is_h5ad:
    if args.output_zarr_path.endswith(".zarr"):
        adata.write_zarr(args.output_zarr_path)
    else:
        adata.write(args.output_zarr_path)
else:
    concatenated_sdata[table_key] = adata
    concatenated_sdata.write(args.output_zarr_path)

    
                   
                   
                   
                   
                   
                   
                   
                   
