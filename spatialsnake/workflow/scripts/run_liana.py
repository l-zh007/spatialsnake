import os
import sys
import inspect
import argparse
import spatialdata as spd
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import sparse
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir in sys.path:
    sys.path.remove(script_dir)
import liana as li
from matplotlib import pyplot as plt
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

plt.rcParams["figure.dpi"] = 50
logger = setup_logger("liana")

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
parser.add_argument("--top_n", type=int, required=False, default=6)
parser.add_argument("--pvalue", type=float, required=False, default=0.05)
parser.add_argument("--source_celltypes", type=str, required=False, default="")
parser.add_argument("--target_celltypes", type=str, required=False, default="")
parser.add_argument("--cell_pairs", type=str, required=False, default="")
parser.add_argument("--pairs", type=str, required=False, default="")
parser.add_argument("--results_csv", type=str, required=False, default="")
parser.add_argument("--top_csv", type=str, required=False, default="")
parser.add_argument("--plot_data_csv", type=str, required=False, default="")
parser.add_argument("--pair_summary_csv", type=str, required=False, default="")
args = parser.parse_args()
output_dir = os.path.dirname(args.output_zarr_path)
if output_dir == "":
    output_dir = "."
os.makedirs(output_dir, exist_ok=True)
results_csv = args.results_csv or os.path.join(output_dir, "liana_results.csv")
top_csv = args.top_csv or os.path.join(output_dir, "liana_top_interactions.csv")
plot_data_csv = args.plot_data_csv or os.path.join(output_dir, "liana_plot_data.csv")
pair_summary_csv = args.pair_summary_csv or os.path.join(output_dir, "liana_pair_summary.csv")
dotplot_path = os.path.join(output_dir, "dotplot.png")
tileplot_path = os.path.join(output_dir, "tileplot.png")
heatmap_path = os.path.join(output_dir, "communication_heatmap.png")
gene_expression_path = os.path.join(output_dir, "gene_expression.png")

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

def split_values(value):
    if value is None:
        return []
    return [x.strip() for x in str(value).replace(";", ",").split(",") if x.strip()]

def parse_pairs(value):
    pairs = []
    for token in split_values(value):
        if "|" in token:
            ligand, receptor = token.split("|", 1)
        elif "->" in token:
            ligand, receptor = token.split("->", 1)
        elif ":" in token:
            ligand, receptor = token.split(":", 1)
        else:
            logger.warning(f"Skip invalid LIANA pair specification: {token}. Use ligand|receptor.")
            continue
        ligand = ligand.strip()
        receptor = receptor.strip()
        if ligand and receptor:
            pairs.append((ligand, receptor))
    return pairs

def parse_cell_pairs(value):
    pairs = []
    for token in split_values(value):
        bidirectional = False
        if "<->" in token:
            source, target = token.split("<->", 1)
            bidirectional = True
        elif "<=>" in token:
            source, target = token.split("<=>", 1)
            bidirectional = True
        elif "->" in token:
            source, target = token.split("->", 1)
        elif "|" in token:
            source, target = token.split("|", 1)
        else:
            logger.warning(f"Skip invalid LIANA cell-pair specification: {token}. Use source|target, source->target, or source<->target.")
            continue
        source = source.strip()
        target = target.strip()
        if source and target:
            pairs.append((source, target))
            if bidirectional:
                pairs.append((target, source))
    unique, seen = [], set()
    for source, target in pairs:
        if (source, target) not in seen:
            unique.append((source, target))
            seen.add((source, target))
    return unique

def is_rank_column(column):
    return column is not None and "rank" in str(column).lower()

def is_pvalue_column(column):
    text = str(column).lower()
    return "pval" in text or "p_value" in text or text.endswith("_p")

def first_existing_column(df, candidates):
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    return None

def save_placeholder_plot(path, message):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.set_axis_off()
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def write_empty_outputs(message):
    empty = pd.DataFrame()
    empty.to_csv(results_csv, index=False)
    empty.to_csv(top_csv, index=False)
    empty.to_csv(plot_data_csv, index=False)
    empty.to_csv(pair_summary_csv, index=False)
    save_placeholder_plot(dotplot_path, message)
    save_placeholder_plot(tileplot_path, message)
    save_placeholder_plot(heatmap_path, message)
    save_placeholder_plot(gene_expression_path, message)

def validate_adata(adata):
    if args.celltype not in adata.obs.columns:
        raise ValueError(
            f"celltype_col '{args.celltype}' was not found in adata.obs. "
            "Please set celltype_col to a valid annotation column."
        )
    labels = adata.obs[args.celltype].astype(str).str.strip()
    labels = labels.replace({"": np.nan, "nan": np.nan, "None": np.nan})
    adata.obs[args.celltype] = pd.Categorical(labels)
    valid_labels = pd.Series(adata.obs[args.celltype]).dropna()
    n_celltypes = valid_labels.nunique()
    if n_celltypes < 2:
        raise ValueError(
            f"LIANA requires at least two valid cell types in '{args.celltype}', "
            f"but found {n_celltypes}."
        )
    counts = valid_labels.value_counts()
    rare = counts[counts < int(args.min_cells)]
    if rare.shape[0] > 0:
        logger.warning(
            "Some cell types have fewer cells than liana_min_cells and may be filtered by LIANA: "
            + ", ".join([f"{idx}({val})" for idx, val in rare.items()])
        )

def sort_liana_results(res_df, magnitude, specificity):
    sort_cols = []
    ascending = []
    if specificity and specificity in res_df.columns:
        sort_cols.append(specificity)
        ascending.append(True if (is_pvalue_column(specificity) or is_rank_column(specificity)) else False)
    if magnitude and magnitude in res_df.columns:
        sort_cols.append(magnitude)
        ascending.append(True if is_rank_column(magnitude) else False)
    if not sort_cols:
        return res_df.copy()
    sorted_df = res_df.copy()
    for column in sort_cols:
        sorted_df[column] = pd.to_numeric(sorted_df[column], errors="coerce")
    sorted_df = sorted_df.sort_values(sort_cols, ascending=ascending, na_position="last")
    return sorted_df

def filter_liana_results(res_df):
    cell_pair_values = parse_cell_pairs(args.cell_pairs)
    source_values = split_values(args.source_celltypes)
    target_values = split_values(args.target_celltypes)
    pair_values = parse_pairs(args.pairs)
    if not cell_pair_values and not source_values and not target_values and not pair_values:
        return res_df.copy(), False

    filtered = res_df.copy()
    requested_filters = []
    cell_pair_filter_applied = False
    used_fallback = False

    if cell_pair_values:
        requested_filters.append("directed cell pair")
        if "source" in filtered.columns and "target" in filtered.columns:
            candidate = filtered.copy()
            candidate["source"] = candidate["source"].astype(str).str.strip()
            candidate["target"] = candidate["target"].astype(str).str.strip()
            cell_pair_set = set(cell_pair_values)
            candidate = candidate[
                candidate.apply(lambda row: (row["source"], row["target"]) in cell_pair_set, axis=1)
            ]
            if candidate.shape[0] > 0:
                filtered = candidate
                cell_pair_filter_applied = True
                logger.info(f"Using {filtered.shape[0]} LIANA rows after directed cell-pair filtering.")
            else:
                logger.warning("None of the requested liana_cell_pairs were found; trying source/target filters or automatic selection.")
                used_fallback = True
        else:
            logger.warning("User requested liana_cell_pairs, but LIANA results have no source/target columns.")
            used_fallback = True

    if not cell_pair_filter_applied and source_values:
        requested_filters.append("source")
        if "source" in filtered.columns:
            filtered["source"] = filtered["source"].astype(str).str.strip()
            filtered = filtered[filtered["source"].isin(source_values)]
        else:
            logger.warning("User requested liana_source_celltypes, but LIANA results have no 'source' column.")
    if not cell_pair_filter_applied and target_values:
        requested_filters.append("target")
        if "target" in filtered.columns:
            filtered["target"] = filtered["target"].astype(str).str.strip()
            filtered = filtered[filtered["target"].isin(target_values)]
        else:
            logger.warning("User requested liana_target_celltypes, but LIANA results have no 'target' column.")
    if pair_values:
        requested_filters.append("ligand-receptor pair")
        before_pair_filter = filtered.copy()
        ligand_col = first_existing_column(filtered, ["ligand", "ligand_complex", "source_genesymbol", "ligand_symbol"])
        receptor_col = first_existing_column(filtered, ["receptor", "receptor_complex", "target_genesymbol", "receptor_symbol"])
        if ligand_col and receptor_col:
            filtered[ligand_col] = filtered[ligand_col].astype(str).str.strip()
            filtered[receptor_col] = filtered[receptor_col].astype(str).str.strip()
            pair_set = set(pair_values)
            filtered = filtered[
                filtered.apply(lambda row: (row[ligand_col], row[receptor_col]) in pair_set, axis=1)
            ]
            if filtered.shape[0] == 0 and before_pair_filter.shape[0] > 0:
                logger.warning("None of the requested liana_pairs were found; using automatic ligand-receptor ranking within the selected cell scope.")
                filtered = before_pair_filter
        else:
            logger.warning("User requested liana_pairs, but ligand/receptor columns were not found in LIANA results.")
    if filtered.shape[0] == 0:
        logger.warning(
            "No LIANA results matched user-defined "
            + ", ".join(requested_filters)
            + "; falling back to automatically selected top interactions."
        )
        return res_df.copy(), True
    logger.info(f"Using {filtered.shape[0]} LIANA rows after user-defined filtering.")
    return filtered, used_fallback

def summarize_cell_pairs(res_df, magnitude, specificity):
    columns = [
        "source", "target", "interaction_count", "significant_interaction_count",
        "mean_score", "best_score", "summary_metric"
    ]
    if res_df is None or res_df.empty or "source" not in res_df.columns or "target" not in res_df.columns:
        return pd.DataFrame(columns=columns)
    summary_input = res_df.copy()
    summary_input["source"] = summary_input["source"].astype(str).str.strip()
    summary_input["target"] = summary_input["target"].astype(str).str.strip()
    if magnitude and magnitude in summary_input.columns:
        summary_input["_score"] = pd.to_numeric(summary_input[magnitude], errors="coerce")
    else:
        summary_input["_score"] = np.nan
    if specificity and specificity in summary_input.columns and is_pvalue_column(specificity):
        pvalues = pd.to_numeric(summary_input[specificity], errors="coerce")
        summary_input["_is_significant"] = pvalues <= float(args.pvalue)
    else:
        summary_input["_is_significant"] = False
    grouped = summary_input.groupby(["source", "target"], observed=True)
    rows = []
    rank_like = is_rank_column(magnitude)
    has_pvalue = specificity and specificity in summary_input.columns and is_pvalue_column(specificity)
    for (source, target), group in grouped:
        scores = group["_score"].dropna()
        n_interactions = int(group.shape[0])
        n_significant = int(group["_is_significant"].sum()) if has_pvalue else np.nan
        mean_score = float(scores.mean()) if not scores.empty else np.nan
        best_score = float(scores.min() if rank_like else scores.max()) if not scores.empty else np.nan
        summary_metric = n_significant if has_pvalue and n_significant > 0 else (
            n_interactions if scores.empty else mean_score
        )
        rows.append({
            "source": source,
            "target": target,
            "interaction_count": n_interactions,
            "significant_interaction_count": n_significant,
            "mean_score": mean_score,
            "best_score": best_score,
            "summary_metric": summary_metric,
        })
    return pd.DataFrame(rows, columns=columns)

def plot_pair_heatmap(pair_summary, path):
    if pair_summary is None or pair_summary.empty:
        save_placeholder_plot(path, "No LIANA source-target summary available")
        return
    metric = (
        "significant_interaction_count"
        if pair_summary["significant_interaction_count"].notna().any()
        and pair_summary["significant_interaction_count"].sum(skipna=True) > 0
        else "summary_metric"
    )
    matrix = pair_summary.pivot_table(
        index="source",
        columns="target",
        values=metric,
        aggfunc="sum",
        fill_value=0,
    )
    if matrix.empty:
        save_placeholder_plot(path, "No LIANA source-target summary available")
        return
    fig_width = min(16, max(5, 0.55 * matrix.shape[1] + 2.5))
    fig_height = min(16, max(4.5, 0.45 * matrix.shape[0] + 2.5))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        matrix,
        cmap=sns.blend_palette(["#edf2f5", "#8fb5c9", "#8b3f5f"], as_cmap=True),
        linewidths=0.5,
        linecolor="white",
        square=False,
        cbar_kws={"label": "Significant interactions" if metric == "n_significant" else "Summary score"},
        ax=ax,
    )
    ax.set_title("LIANA source-target communication summary", loc="left", fontweight="bold")
    ax.set_xlabel("Receiver cell type")
    ax.set_ylabel("Sender cell type")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def extract_interaction_genes(res_df, limit=12):
    genes = []
    for column in ["ligand", "receptor", "ligand_complex", "receptor_complex"]:
        if column not in res_df.columns:
            continue
        for value in res_df[column].dropna().astype(str):
            tokens = [value]
            if "_" in value:
                tokens.extend([part for part in value.split("_") if part])
            for token in tokens:
                token = token.strip()
                if token and token not in genes:
                    genes.append(token)
                if len(genes) >= limit:
                    return genes
    return genes

def expression_source(adata):
    use_raw = parse_bool(args.use_raw) and adata.raw is not None
    if use_raw:
        return adata.raw, pd.Index(adata.raw.var_names.astype(str))
    return adata, pd.Index(adata.var_names.astype(str))

def plot_gene_expression(adata, genes, celltypes, path):
    expr_obj, var_names = expression_source(adata)
    genes = [gene for gene in genes if gene in var_names]
    celltypes = [cell_type for cell_type in celltypes if cell_type in set(adata.obs[args.celltype].astype(str))]
    if not genes or not celltypes:
        save_placeholder_plot(path, "No selected ligand-receptor genes were found in the expression matrix")
        return
    means, fractions = [], []
    obs_labels = adata.obs[args.celltype].astype(str)
    for cell_type in celltypes:
        mask = obs_labels.eq(cell_type).to_numpy()
        matrix = expr_obj[mask, genes].X
        if sparse.issparse(matrix):
            means.append(np.asarray(matrix.mean(axis=0)).ravel())
            fractions.append(np.asarray((matrix > 0).mean(axis=0)).ravel())
        else:
            matrix = np.asarray(matrix)
            means.append(matrix.mean(axis=0))
            fractions.append((matrix > 0).mean(axis=0))
    means = np.asarray(means, dtype=float)
    fractions = np.asarray(fractions, dtype=float)
    minimum = np.nanmin(means, axis=0, keepdims=True)
    span = np.nanmax(means, axis=0, keepdims=True) - minimum
    scaled = np.divide(means - minimum, span, out=np.zeros_like(means), where=span > 0)
    width = min(18, max(6, 0.45 * len(celltypes) + 2.5))
    height = min(12, max(4, 0.34 * len(genes) + 2.0))
    fig, ax = plt.subplots(figsize=(width, height))
    x, y = np.meshgrid(np.arange(len(celltypes)), np.arange(len(genes)))
    scatter = ax.scatter(
        x.ravel(),
        y.ravel(),
        s=(20 + 110 * fractions.T).ravel(),
        c=scaled.T.ravel(),
        cmap="viridis",
        vmin=0,
        vmax=1,
        edgecolors="white",
        linewidths=0.35,
    )
    ax.set_xticks(range(len(celltypes)), celltypes, rotation=45, ha="right")
    ax.set_yticks(range(len(genes)), genes)
    ax.set_xlabel("Cell type")
    ax.set_ylabel("")
    ax.set_title("Selected ligand-receptor gene expression", loc="left", fontweight="bold")
    ax.invert_yaxis()
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.035)
    colorbar.set_label("Mean expression (scaled per gene)")
    handles = [
        ax.scatter([], [], s=20 + 110 * fraction, color="#7a9bb8", edgecolors="white", label=f"{int(fraction * 100)}%")
        for fraction in (0.25, 0.5, 0.75)
    ]
    ax.legend(handles=handles, title="Cells expressed", bbox_to_anchor=(1.12, 1), loc="upper left")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def run_method(adata, method_name, uns_key):
    method = getattr(li.mt, method_name, None)
    if method is None:
        available = [m for m in dir(li.mt) if not m.startswith("_")]
        raise ValueError(f"Unsupported method: {method_name}. Available methods: {available}")
    use_raw = parse_bool(args.use_raw)
    if use_raw and adata.raw is None:
        logger.warning("adata.raw is not initialized; switching use_raw to False")
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
    if not isinstance(res_df, pd.DataFrame) or res_df.shape[0] == 0:
        logger.warning("Skip LIANA plotting/export: LIANA returned no tabular results.")
        write_empty_outputs("LIANA returned no tabular ligand-receptor results")
        return
    magnitude, specificity = get_method_config(method_name, res_df)
    if magnitude is None:
        logger.warning("Skip LIANA plotting: no compatible LIANA score column was found.")
        res_df.to_csv(results_csv, index=False)
        res_df.to_csv(top_csv, index=False)
        res_df.to_csv(plot_data_csv, index=False)
        pd.DataFrame().to_csv(pair_summary_csv, index=False)
        save_placeholder_plot(dotplot_path, "No compatible LIANA score column was found")
        save_placeholder_plot(tileplot_path, "No compatible LIANA score column was found")
        save_placeholder_plot(heatmap_path, "No compatible LIANA score column was found")
        save_placeholder_plot(gene_expression_path, "No compatible LIANA score column was found")
        return
    adata.obs[args.celltype] = pd.Categorical(adata.obs[args.celltype].astype(str).str.strip())
    celltypes = adata.obs[args.celltype].cat.categories.astype(str).tolist()
    if len(celltypes) == 0:
        return
    plot_uns_key = uns_key

    res_df.to_csv(results_csv, index=False)
    if "source" in res_df.columns and "target" in res_df.columns:
        valid_cells = set(celltypes)
        res_df = res_df.copy()
        res_df["source"] = res_df["source"].astype(str).str.strip()
        res_df["target"] = res_df["target"].astype(str).str.strip()
        res_df = res_df[res_df["source"].isin(valid_cells) & res_df["target"].isin(valid_cells)]
        if res_df.shape[0] == 0:
            logger.warning("Skip LIANA plotting: no source/target pairs matched current celltype categories")
            pd.DataFrame().to_csv(top_csv, index=False)
            pd.DataFrame().to_csv(plot_data_csv, index=False)
            pd.DataFrame().to_csv(pair_summary_csv, index=False)
            save_placeholder_plot(dotplot_path, "No LIANA source/target pairs matched current cell type categories")
            save_placeholder_plot(tileplot_path, "No LIANA source/target pairs matched current cell type categories")
            save_placeholder_plot(heatmap_path, "No LIANA source/target pairs matched current cell type categories")
            save_placeholder_plot(gene_expression_path, "No LIANA source/target pairs matched current cell type categories")
            return

    sorted_df = sort_liana_results(res_df, magnitude, specificity)
    top_n = max(1, int(args.top_n))
    top_limit = min(max(top_n * max(1, len(celltypes)), top_n), len(sorted_df))
    top_df = sorted_df.head(top_limit).copy()
    top_df.to_csv(top_csv, index=False)

    filtered_df, used_fallback = filter_liana_results(sorted_df)
    filtered_df = sort_liana_results(filtered_df, magnitude, specificity)
    pair_summary = summarize_cell_pairs(filtered_df, magnitude, specificity)
    pair_summary.to_csv(pair_summary_csv, index=False)
    plot_pair_heatmap(pair_summary, heatmap_path)
    pvalue_filter_applied = False
    if specificity and specificity in filtered_df.columns and is_pvalue_column(specificity):
        sig_df = filtered_df[pd.to_numeric(filtered_df[specificity], errors="coerce") <= float(args.pvalue)].copy()
        if sig_df.shape[0] > 0:
            filtered_df = sig_df
            pvalue_filter_applied = True
        else:
            logger.warning(
                f"No LIANA interactions passed liana_pvalue <= {args.pvalue}; "
                "plots will show the best available non-significant interactions."
            )
    plot_limit = min(max(top_n * max(1, len(celltypes)), top_n), len(filtered_df))
    res_df = filtered_df.head(plot_limit).copy()
    res_df.to_csv(plot_data_csv, index=False)

    if res_df.shape[0] == 0:
        logger.warning("Skip LIANA plotting: no interactions available after filtering.")
        save_placeholder_plot(dotplot_path, "No LIANA interactions available after filtering")
        save_placeholder_plot(tileplot_path, "No LIANA interactions available after filtering")
        save_placeholder_plot(gene_expression_path, "No LIANA interactions available after filtering")
        return
    if used_fallback:
        logger.warning("LIANA plots use automatic top interactions because custom filters produced no matches.")

    plot_uns_key = f"{uns_key}__plot"
    adata.uns[plot_uns_key] = res_df
    source_labels = celltypes
    target_labels = celltypes
    if "source" in res_df.columns:
        available_sources = set(res_df["source"].astype(str).unique().tolist())
        source_labels = [s for s in source_labels if s in available_sources]
    if "target" in res_df.columns:
        available_targets = set(res_df["target"].astype(str).unique().tolist())
        target_labels = [t for t in target_labels if t in available_targets]
    if len(source_labels) == 0 or len(target_labels) == 0:
        logger.warning("Skip LIANA plotting: no source/target labels remain after filtering.")
        save_placeholder_plot(dotplot_path, "No source/target labels remain after LIANA filtering")
        save_placeholder_plot(tileplot_path, "No source/target labels remain after LIANA filtering")
        save_placeholder_plot(gene_expression_path, "No source/target labels remain after LIANA filtering")
        return
    filter_fun = None
    inverse_size = False
    orderby_ascending = False
    inverse_score = False
    if specificity and specificity in res_df.columns and is_pvalue_column(specificity) and pvalue_filter_applied:
        def filter_fun(x):
            return x[specificity] <= float(args.pvalue)
        inverse_size = True
        orderby_ascending = True
        if method_name == "cellphonedb":
            inverse_score = True
    elif specificity and specificity in res_df.columns and is_rank_column(specificity):
        orderby_ascending = True
    size_key = specificity if specificity else magnitude
    plot_celltypes = [cell_type for cell_type in celltypes if cell_type in set(source_labels + target_labels)]
    dotplot_fig = li.pl.dotplot(
        adata=adata,
        colour=magnitude,
        size=size_key,
        inverse_size=inverse_size,
        source_labels=source_labels,
        target_labels=target_labels,
        figure_size=(max(12, len(target_labels) * 1.2), max(8, len(source_labels) * 0.8)),
        filter_fun=filter_fun,
        uns_key=plot_uns_key
    )
    dotplot_fig.save(dotplot_path, dpi=300)
    if {"ligand_means", "receptor_means", "ligand_props", "receptor_props"}.issubset(res_df.columns):
        tile_fill = "means"
        tile_label = "props"
        tile_label_fun = lambda x: f"{x:.2f}"
        tile_orderby = size_key
        tile_orderby_ascending = orderby_ascending
    else:
        tile_fill = magnitude
        tile_label = size_key
        tile_label_fun = lambda x: f"{x:.2f}"
        tile_orderby = size_key
        tile_orderby_ascending = orderby_ascending
    my_plot = li.pl.tileplot(
        adata=adata,
        fill=tile_fill,
        label=tile_label,
        label_fun=tile_label_fun,
        top_n=max(1, int(args.top_n)),
        orderby=tile_orderby,
        orderby_ascending=tile_orderby_ascending,
        source_labels=source_labels,
        target_labels=target_labels,
        uns_key=plot_uns_key,
        source_title="Ligand",
        target_title="Receptor",
        figure_size=(max(10, len(target_labels) * 1.1), max(8, len(source_labels) * 0.8))
    )
    my_plot.save(tileplot_path, dpi=300)
    plot_gene_expression(
        adata,
        extract_interaction_genes(res_df, limit=12),
        plot_celltypes,
        gene_expression_path,
    )
    try:
        li.pl.circle_plot(
            adata,
            groupby=args.celltype,
            score_key=magnitude,
            inverse_score=inverse_score,
            source_labels=source_labels,
            target_labels=target_labels,
            filter_fun=filter_fun,
            top_n=max(1, int(args.top_n)),
            orderby=size_key,
            orderby_ascending=orderby_ascending,
            pivot_mode="counts",
            figure_size=(10, 10),
            uns_key=plot_uns_key
        )
    except Exception as first_error:
        try:
            li.pl.circle_plot(
                adata,
                groupby=args.celltype,
                score_key=magnitude,
                inverse_score=inverse_score,
                source_labels=celltypes,
                target_labels=celltypes,
                filter_fun=None,
                top_n=max(1, int(args.top_n)),
                orderby=size_key,
                orderby_ascending=orderby_ascending,
                pivot_mode="counts",
                figure_size=(10, 10),
                uns_key=plot_uns_key
            )
        except Exception as second_error:
            logger.warning(
                "Skip circle plot because LIANA could not construct a stable network "
                f"from the selected results: {first_error}; fallback error: {second_error}"
            )
            plt.close()
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
    log_step(logger, 1, 4, "loading AnnData input")
    adata = sc.read_h5ad(args.input_dir)
else:
    log_step(logger, 1, 4, "loading SpatialData input")
    concatenated_sdata = spd.read_zarr(args.input_dir)
    table_key = next(iter(concatenated_sdata.tables.keys()))
    adata = concatenated_sdata[table_key]
logger.info(f"Loaded {adata.n_obs} observations and {adata.n_vars} genes")
validate_adata(adata)

log_step(logger, 2, 4, f"running LIANA method={method_name}")
run_method(adata, method_name, uns_key)
log_step(logger, 3, 4, "saving LIANA plots")
plot_results(adata, method_name, uns_key, output_dir)

log_step(logger, 4, 4, f"saving LIANA output to {args.output_zarr_path}")
if input_is_h5ad:
    if args.output_zarr_path.endswith(".zarr"):
        adata.write_zarr(args.output_zarr_path)
    else:
        adata.write(args.output_zarr_path)
else:
    concatenated_sdata[table_key] = adata
    concatenated_sdata.write(args.output_zarr_path)
logger.info("LIANA module completed")

    
                   
                   
                   
                   
                   
                   
                   
                   
