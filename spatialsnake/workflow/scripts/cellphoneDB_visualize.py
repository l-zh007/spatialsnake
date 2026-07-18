import argparse
import os
import re

import anndata as ad
import ktplotspy as kpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from plotnine import element_text, guides, theme
from scipy import sparse

from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("cellphonedb_visualize")
parser = argparse.ArgumentParser(description="Visualize CellPhoneDB results")
parser.add_argument("--input_dir", required=True, help="CellPhoneDB input h5ad")
parser.add_argument("--output_zarr_path", required=True, help="Output heatmap path")
parser.add_argument("--sample_id", required=True, help="Sample identifier")
parser.add_argument("--type", default="", help="Input technology (kept for workflow compatibility)")
parser.add_argument("--output_name", default="", help="Suffix used by CellPhoneDB result files")
parser.add_argument("--cell_pairs", default="", help="Comma-separated directed or bidirectional cell pairs, e.g. A|B,A<->C")
parser.add_argument("--cell_type1", default="", help="Optional sender cell type")
parser.add_argument("--cell_type2", default="", help="Optional receiver cell type")
parser.add_argument("--gene_family", default="", help="Optional ktplotspy gene family")
parser.add_argument("--cpdb_pathway", default="", help="Optional CellPhoneDB interaction classification")
parser.add_argument("--interaction_pairs", default="", help="Comma-separated ligand-receptor pairs")
parser.add_argument("--cpdb_genes", default="", help="Comma-separated genes for the expression plot")
parser.add_argument("--celltype", default="", help="Cell-type annotation column")
parser.add_argument("--method", default="", choices=["", "statistical", "degs"], help="CellPhoneDB method")
parser.add_argument("--pvalue", type=float, default=0.05, help="Significance threshold")
args = parser.parse_args()


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


def split_values(value):
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def split_pair_values(value):
    return [item.strip() for item in re.split(r"[,;]", str(value or "")) if item.strip()]


def safe_label(value):
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.")
    return label or "selected"


def output_paths(output_dir, sample_id, cell_type1, cell_type2, suffix, legacy_name=None, include_pair_specific=False):
    paths = []
    if legacy_name:
        paths.append(os.path.join(output_dir, legacy_name))
    if include_pair_specific:
        pair_label = f"{safe_label(cell_type1)}_to_{safe_label(cell_type2)}"
        paths.append(os.path.join(output_dir, f"{sample_id}_{pair_label}_{suffix}"))
    seen, unique = set(), []
    for path in paths:
        if path not in seen:
            unique.append(path)
            seen.add(path)
    return unique


def pick_file(base_dir, candidates, required=True):
    for candidate in candidates:
        path = os.path.join(base_dir, candidate)
        if os.path.isfile(path):
            return path
    if required:
        raise FileNotFoundError(f"Missing CellPhoneDB result file; tried: {candidates}")
    return None


def resolve_celltype_key(adata, requested):
    if requested and requested in adata.obs.columns:
        return requested
    stored = adata.uns.get("cellphonedb", {}).get("celltype_col", "")
    if stored in adata.obs.columns:
        return stored
    for candidate in ["celltype", "cell_type", "celltypes", "cell_type_annotation"]:
        if candidate in adata.obs.columns:
            return candidate
    raise ValueError("No cell-type annotation column was found in adata.obs")


def pair_columns(frame):
    return [column for column in frame.columns if "|" in str(column)]


def pair_statistics(significance, degs_analysis, alpha):
    columns = pair_columns(significance)
    if not columns:
        raise ValueError("No directed cell-type pair columns were found in CellPhoneDB results")
    numeric = significance[columns].apply(pd.to_numeric, errors="coerce")
    significant = numeric.gt(0) if degs_analysis else numeric.lt(alpha)
    counts = significant.sum(axis=0).sort_values(ascending=False)
    return columns, significant, counts


def parse_pair_token(token):
    token = str(token or "").strip()
    for separator, bidirectional in [("<->", True), ("<=>", True), ("->", False), ("|", False)]:
        if separator in token:
            left, right = [part.strip() for part in token.split(separator, 1)]
            if not left or not right:
                return []
            pairs = [(left, right)]
            if bidirectional:
                pairs.append((right, left))
            return pairs
    return []


def pair_names(cell_type1, cell_type2):
    return f"{cell_type1}|{cell_type2}"


def select_cell_pairs(significance, requested_pairs, requested_1, requested_2, degs_analysis, alpha):
    columns, significant, counts = pair_statistics(significance, degs_analysis, alpha)
    available = set(columns)
    selected = []
    selection_source = "automatic"

    requested_pair_tokens = split_pair_values(requested_pairs)
    if requested_pair_tokens:
        selection_source = "cell_pairs"
        for token in requested_pair_tokens:
            parsed_pairs = parse_pair_token(token)
            if not parsed_pairs:
                logger.warning("Could not parse requested cell pair '%s'; use A|B, A->B, or A<->B", token)
                continue
            for cell_type1, cell_type2 in parsed_pairs:
                pair = pair_names(cell_type1, cell_type2)
                if pair in available:
                    selected.append((pair, cell_type1, cell_type2, selection_source))
                else:
                    logger.warning("Requested cell pair '%s' was not found in CellPhoneDB results", pair)

    requested_1, requested_2 = str(requested_1 or "").strip(), str(requested_2 or "").strip()
    if not selected and (requested_1 or requested_2):
        selection_source = "cell_type1_cell_type2"
        exact = pair_names(requested_1, requested_2) if requested_1 and requested_2 else ""
        if exact in available:
            selected.append((exact, requested_1, requested_2, selection_source))
        else:
            candidates = columns
            if requested_1:
                candidates = [column for column in candidates if column.split("|", 1)[0] == requested_1]
            if requested_2:
                candidates = [column for column in candidates if column.split("|", 1)[1] == requested_2]
            if candidates:
                pair = counts[counts.index.isin(candidates)].index[0]
                cell_type1, cell_type2 = pair.split("|", 1)
                logger.warning(
                    "Exact requested cell pair '%s' was not found; plotting '%s' instead",
                    exact or f"{requested_1 or '*'}|{requested_2 or '*'}",
                    pair,
                )
                selected.append((pair, cell_type1, cell_type2, selection_source))
            else:
                logger.warning(
                    "Requested cell type filter '%s' did not match any directed pair; using automatic selection",
                    exact or f"{requested_1 or '*'}|{requested_2 or '*'}",
                )

    if not selected:
        pair = counts.index[0]
        cell_type1, cell_type2 = pair.split("|", 1)
        selected.append((pair, cell_type1, cell_type2, "automatic"))

    unique, seen = [], set()
    for pair, cell_type1, cell_type2, source in selected:
        if pair in seen:
            continue
        seen.add(pair)
        unique.append(
            {
                "pair": pair,
                "cell_type1": cell_type1,
                "cell_type2": cell_type2,
                "source": source,
                "significant_rows": significant[pair],
                "n_significant": int(counts[pair]),
            }
        )
        logger.info("Selected directed cell pair: %s (%d significant interactions)", pair, int(counts[pair]))
    return unique


def align_column(source, target, column, default=np.nan):
    if source is None or column not in source.columns:
        return pd.Series(default, index=target.index, dtype=float)
    key = "id_cp_interaction"
    if key in source.columns and key in target.columns:
        values = pd.to_numeric(source.set_index(key)[column], errors="coerce")
        return target[key].map(values).astype(float)
    values = pd.to_numeric(source[column], errors="coerce").reset_index(drop=True)
    return values.reindex(range(len(target))).set_axis(target.index)


def rank_interactions(means, significance, scores, pair, significant_rows, degs_analysis):
    ranked = means.loc[significant_rows.to_numpy()].copy()
    if ranked.empty:
        return ranked
    ranked["_mean"] = align_column(means, ranked, pair)
    ranked["_score"] = align_column(scores, ranked, pair)
    ranked["_significance"] = align_column(significance, ranked, pair)
    ranked["_pvalue_sort"] = 0.0 if degs_analysis else ranked["_significance"].fillna(1.0)
    ranked["_score_sort"] = ranked["_score"].fillna(-1.0)
    ranked["_mean_sort"] = ranked["_mean"].fillna(-1.0)
    return ranked.sort_values(
        ["_score_sort", "_mean_sort", "_pvalue_sort"],
        ascending=[False, False, True],
    )


def interacting_pairs_from_table(table, limit):
    if "interacting_pair" not in table.columns:
        return []
    return table["interacting_pair"].dropna().astype(str).drop_duplicates().head(limit).tolist()


def filter_requested_pairs(ranked, requested):
    requested_pairs = split_values(requested)
    if not requested_pairs or "interacting_pair" not in ranked.columns:
        return ranked
    selected = ranked[ranked["interacting_pair"].astype(str).isin(requested_pairs)]
    if selected.empty:
        logger.warning("None of the requested interaction_pairs were significant for the selected cell pair; using automatic ranking")
        return ranked
    return selected


def select_focused_interactions(ranked, requested_pairs, requested_pathway):
    if requested_pairs:
        return filter_requested_pairs(ranked, requested_pairs), "user-selected interactions"
    if "classification" not in ranked.columns:
        return ranked, "top interactions"

    classifications = ranked["classification"].fillna("").astype(str).str.strip()
    valid = ranked[classifications.ne("") & classifications.str.lower().ne("nan")]
    if valid.empty:
        return ranked, "top interactions"

    if requested_pathway:
        mask = valid["classification"].astype(str).str.contains(str(requested_pathway), case=False, regex=False)
        selected = valid[mask]
        if not selected.empty:
            return selected, str(selected["classification"].iloc[0])
        logger.warning("cpdb_pathway '%s' was not found; using the most represented significant pathway", requested_pathway)

    pathway = valid["classification"].value_counts().index[0]
    return valid[valid["classification"] == pathway], str(pathway)


def save_plotnine(plot, paths, width, height):
    if isinstance(paths, str):
        paths = [paths]
    plot = plot + guides(size="none", stroke="none", alpha="none") + theme(
        legend_box="vertical",
        legend_position="right",
        axis_text_x=element_text(angle=45, hjust=1),
    )
    for path in paths:
        plot.save(path, dpi=300, width=width, height=height, units="in", limitsize=False, verbose=False)
    plt.close("all")


def plot_gene_expression(adata, celltype_key, genes, output_paths):
    genes = [gene for gene in genes if gene in adata.var_names]
    if not genes:
        logger.info("Skipping gene-expression plot because no selected genes were found in adata.var_names")
        return

    cell_types = adata.obs[celltype_key].dropna().astype(str).drop_duplicates().tolist()
    means, fractions = [], []
    for cell_type in cell_types:
        matrix = adata[adata.obs[celltype_key].astype(str).eq(cell_type), genes].X
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

    width = min(18, max(6, 0.42 * len(cell_types) + 2.5))
    height = min(12, max(4, 0.34 * len(genes) + 2.0))
    fig, ax = plt.subplots(figsize=(width, height))
    x, y = np.meshgrid(np.arange(len(cell_types)), np.arange(len(genes)))
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
    ax.set_xticks(range(len(cell_types)), cell_types, rotation=45, ha="right")
    ax.set_yticks(range(len(genes)), genes)
    ax.set_xlabel("Cell type")
    ax.set_ylabel("")
    ax.set_title("Selected ligand-receptor gene expression", loc="left", fontweight="bold")
    ax.invert_yaxis()
    ax.grid(False)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.035)
    colorbar.set_label("Mean expression (scaled per gene)")
    handles = [
        ax.scatter([], [], s=20 + 110 * fraction, color="#7a9bb8", edgecolors="white", label=f"{int(fraction * 100)}%")
        for fraction in (0.25, 0.5, 0.75)
    ]
    ax.legend(handles=handles, title="Cells expressed", bbox_to_anchor=(1.12, 1), loc="upper left")
    if isinstance(output_paths, str):
        output_paths = [output_paths]
    for output_path in output_paths:
        fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def export_interaction_rows(table, pair_info, panel, limit, selected_genes, degs_analysis):
    rows = []
    if table is None or table.empty:
        rows.append(
            {
                "selected_pair": pair_info["pair"],
                "sender": pair_info["cell_type1"],
                "receiver": pair_info["cell_type2"],
                "selection_source": pair_info["source"],
                "plot_panel": panel,
                "rank": "",
                "n_significant_for_pair": pair_info["n_significant"],
                "interacting_pair": "",
                "classification": "",
                "gene_a": "",
                "gene_b": "",
                "mean": "",
                "interaction_score": "",
                "pvalue": "",
                "relevant_interaction": "",
                "selected_genes": ";".join(selected_genes),
            }
        )
        return rows

    for rank, (_, row) in enumerate(table.head(limit).iterrows(), start=1):
        significance_value = row.get("_significance", "")
        rows.append(
            {
                "selected_pair": pair_info["pair"],
                "sender": pair_info["cell_type1"],
                "receiver": pair_info["cell_type2"],
                "selection_source": pair_info["source"],
                "plot_panel": panel,
                "rank": rank,
                "n_significant_for_pair": pair_info["n_significant"],
                "interacting_pair": row.get("interacting_pair", ""),
                "classification": row.get("classification", ""),
                "gene_a": row.get("gene_a", ""),
                "gene_b": row.get("gene_b", ""),
                "mean": row.get("_mean", ""),
                "interaction_score": row.get("_score", ""),
                "pvalue": "" if degs_analysis else significance_value,
                "relevant_interaction": significance_value if degs_analysis else "",
                "selected_genes": ";".join(selected_genes),
            }
        )
    return rows


def plot_microenvironments(path, output_path):
    if not path or not os.path.isfile(path):
        return
    microenvs = pd.read_csv(path, sep=None, engine="python", dtype=str)
    if microenvs.shape[1] < 2:
        return
    microenvs = microenvs.iloc[:, :2].dropna().drop_duplicates()
    microenvs.columns = ["cell_type", "microenvironment"]
    matrix = pd.crosstab(microenvs["cell_type"], microenvs["microenvironment"]).gt(0).astype(int)
    if matrix.empty:
        return
    fig, ax = plt.subplots(
        figsize=(min(14, max(5, 0.55 * matrix.shape[1] + 2)), min(14, max(4, 0.35 * matrix.shape[0] + 2)))
    )
    sns.heatmap(
        matrix,
        cmap=sns.color_palette(["#f1f3f5", "#315d7d"], as_cmap=True),
        cbar=False,
        linewidths=0.5,
        linecolor="white",
        square=False,
        ax=ax,
    )
    ax.set_xlabel("Microenvironment")
    ax.set_ylabel("Cell type")
    ax.set_title("Cell type–microenvironment membership", loc="left", fontweight="bold")
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


log_step(logger, 1, 5, "loading CellPhoneDB input and result files")
adata = ad.read_h5ad(args.input_dir)
run_metadata = adata.uns.get("cellphonedb", {})
output_dir = os.path.dirname(args.output_zarr_path)
base_dir = os.path.join(output_dir, "cellphonedb_output")
method = args.method or str(run_metadata.get("method", "statistical"))
degs_analysis = method == "degs"
prefix = "degs_analysis" if degs_analysis else "statistical_analysis"
output_name = args.output_name.strip() or str(run_metadata.get("output_name", "")).strip() or args.sample_id

means_path = pick_file(base_dir, [f"{prefix}_means_{output_name}.txt"])
if degs_analysis:
    significance_path = pick_file(base_dir, [f"{prefix}_relevant_interactions_{output_name}.txt"])
else:
    significance_path = pick_file(base_dir, [f"{prefix}_pvalues_{output_name}.txt"])
scores_path = pick_file(base_dir, [f"{prefix}_interaction_scores_{output_name}.txt"], required=False)
decon_path = pick_file(
    base_dir,
    [f"{prefix}_deconvoluted_{output_name}.txt", f"{prefix}_deconvoluted_percents_{output_name}.txt"],
    required=False,
)

means = pd.read_csv(means_path, sep="\t")
significance = pd.read_csv(significance_path, sep="\t")
scores = pd.read_csv(scores_path, sep="\t") if scores_path else None
deconvoluted = pd.read_csv(decon_path, sep="\t") if decon_path else None
celltype_key = resolve_celltype_key(adata, args.celltype)
selected_pairs = select_cell_pairs(
    significance,
    args.cell_pairs,
    args.cell_type1,
    args.cell_type2,
    degs_analysis,
    args.pvalue,
)

log_step(logger, 2, 5, "saving the interaction-summary heatmap")
n_cell_types = len(adata.obs[celltype_key].dropna().unique())
heatmap_size = min(14, max(5, 0.45 * n_cell_types + 2.5))
heatmap_tables = kpy.plot_cpdb_heatmap(
    pvals=significance,
    degs_analysis=degs_analysis,
    alpha=args.pvalue,
    return_tables=True,
)
fig, ax = plt.subplots(figsize=(heatmap_size, heatmap_size))
sns.heatmap(
    heatmap_tables["count_network"],
    cmap=sns.blend_palette(["#edf2f5", "#8fb5c9", "#8b3f5f"], as_cmap=True),
    linewidths=0.5,
    linecolor="white",
    square=True,
    cbar_kws={"label": "Significant interactions", "shrink": 0.75},
    ax=ax,
)
ax.set_title(f"{args.sample_id} | Significant interactions", loc="left", fontweight="bold")
ax.set_xlabel("Receiver cell type")
ax.set_ylabel("Sender cell type")
ax.tick_params(axis="x", rotation=45)
for label in ax.get_xticklabels():
    label.set_horizontalalignment("right")
fig.savefig(args.output_zarr_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

log_step(logger, 3, 5, "saving interaction dot plots")
selected_rows = []
write_pair_specific = len(selected_pairs) > 1
allowed_families = {"chemokines", "th1", "th2", "th17", "treg", "costimulatory", "coinhibitory"}
gene_family = args.gene_family.strip().lower()

for pair_index, pair_info in enumerate(selected_pairs, start=1):
    pair = pair_info["pair"]
    cell_type1 = pair_info["cell_type1"]
    cell_type2 = pair_info["cell_type2"]
    significant_rows = pair_info["significant_rows"]
    ranked = rank_interactions(means, significance, scores, pair, significant_rows, degs_analysis)
    ranked = filter_requested_pairs(ranked, args.interaction_pairs)
    top_interactions = interacting_pairs_from_table(ranked, 20)
    is_primary_pair = pair_index == 1
    include_pair_specific = write_pair_specific or not is_primary_pair

    if top_interactions:
        dot_height = min(12, max(4.5, 0.32 * len(top_interactions) + 2.0))
        dot_plot = kpy.plot_cpdb(
            adata=adata,
            cell_type1=cell_type1,
            cell_type2=cell_type2,
            means=means,
            pvals=significance,
            celltype_key=celltype_key,
            interacting_pairs=top_interactions,
            degs_analysis=degs_analysis,
            alpha=args.pvalue,
            standard_scale=False,
            max_size=6,
            highlight_size=1.0,
            title=f"{cell_type1} → {cell_type2}",
        )
        save_plotnine(
            dot_plot,
            output_paths(
                output_dir,
                args.sample_id,
                cell_type1,
                cell_type2,
                "dot_plot.png",
                legacy_name=f"{args.sample_id}_dot_plot.png" if is_primary_pair else None,
                include_pair_specific=include_pair_specific,
            ),
            6.5,
            dot_height,
        )
    else:
        logger.info("Skipping interaction dot plot for %s because no significant interactions were found", pair)

    focused_table, focused_label = select_focused_interactions(ranked, args.interaction_pairs, args.cpdb_pathway)
    focused_pairs = interacting_pairs_from_table(focused_table, 20)
    focused_kwargs = {}
    if not args.interaction_pairs and not args.cpdb_pathway and gene_family in allowed_families:
        focused_kwargs["gene_family"] = gene_family
        focused_label = gene_family
    elif focused_pairs:
        focused_kwargs["interacting_pairs"] = focused_pairs

    if focused_kwargs:
        focused_plot = kpy.plot_cpdb(
            adata=adata,
            cell_type1=cell_type1,
            cell_type2=cell_type2,
            means=means,
            pvals=significance,
            celltype_key=celltype_key,
            degs_analysis=degs_analysis,
            alpha=args.pvalue,
            standard_scale=False,
            max_size=6,
            highlight_size=1.0,
            title=f"Focused signaling | {focused_label}",
            **focused_kwargs,
        )
        focused_count = len(focused_pairs) if focused_pairs else 12
        save_plotnine(
            focused_plot,
            output_paths(
                output_dir,
                args.sample_id,
                cell_type1,
                cell_type2,
                "dot_family_plot.png",
                legacy_name=f"{args.sample_id}_dot_family_plot.png" if is_primary_pair else None,
                include_pair_specific=include_pair_specific,
            ),
            6.5,
            min(11, max(4.5, 0.32 * focused_count + 2.0)),
        )

    selected_genes = split_values(args.cpdb_genes)
    if not selected_genes:
        gene_columns = [column for column in ["gene_a", "gene_b"] if column in ranked.columns]
        if gene_columns:
            selected_genes = (
                ranked.head(20)[gene_columns]
                .stack()
                .dropna()
                .astype(str)
                .loc[lambda values: values.str.strip().ne("")]
                .drop_duplicates()
                .head(12)
                .tolist()
            )

    selected_rows.extend(export_interaction_rows(ranked, pair_info, "dot_plot", 20, selected_genes, degs_analysis))
    selected_rows.extend(export_interaction_rows(focused_table, pair_info, "dot_family_plot", 20, selected_genes, degs_analysis))

    log_step(logger, 4, 5, "saving focused expression and chord summaries")
    plot_gene_expression(
        adata,
        celltype_key,
        selected_genes,
        output_paths(
            output_dir,
            args.sample_id,
            cell_type1,
            cell_type2,
            "gene_expression.png",
            legacy_name=f"{args.sample_id}_gene_expression.png" if is_primary_pair else None,
            include_pair_specific=include_pair_specific,
        ),
    )

    chord_pairs = interacting_pairs_from_table(ranked, 8)
    chord_table = ranked[ranked["interacting_pair"].astype(str).isin(chord_pairs)] if chord_pairs and "interacting_pair" in ranked.columns else ranked.iloc[0:0]
    selected_rows.extend(export_interaction_rows(chord_table, pair_info, "chord_plot", 8, selected_genes, degs_analysis))
    if chord_pairs and deconvoluted is not None:
        try:
            chord = kpy.plot_cpdb_chord(
                adata=adata,
                cell_type1=cell_type1,
                cell_type2=cell_type2,
                means=means,
                pvals=significance,
                deconvoluted=deconvoluted,
                celltype_key=celltype_key,
                interaction=chord_pairs,
                degs_analysis=degs_analysis,
                alpha=args.pvalue,
                link_kwargs={"direction": 1, "allow_twist": True, "r1": 95, "r2": 90},
                sector_text_kwargs={"color": "black", "size": 8, "r": 105, "adjust_rotation": True},
                legend_kwargs={"loc": "center", "bbox_to_anchor": (1, 1), "fontsize": 7},
                link_offset=1,
            )
            for path in output_paths(
                output_dir,
                args.sample_id,
                cell_type1,
                cell_type2,
                "chord_plot.png",
                legacy_name=f"{args.sample_id}_chord_plot.png" if is_primary_pair else None,
                include_pair_specific=include_pair_specific,
            ):
                chord.savefig(path, dpi=300)
            plt.close("all")
        except Exception as error:
            logger.warning("Chord plot for %s was skipped: %s", pair, error)

pd.DataFrame(selected_rows).to_csv(os.path.join(output_dir, f"{args.sample_id}_selected_interactions.csv"), index=False)

microenvs_file_path = str(run_metadata.get("microenvs_file_path", ""))
plot_microenvironments(
    microenvs_file_path,
    os.path.join(output_dir, f"{args.sample_id}_microenvironment.png"),
)

log_step(logger, 5, 5, "CellPhoneDB visualization completed")
