import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as spd
from matplotlib.lines import Line2D

from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.logging_utils import log_step, setup_logger
from spatialsnake.workflow.function.plot import plot_cell2location_dotplot


logger = setup_logger("RCTD_visualize")


def _natural_key(value):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]


def _validate_pdf_path(path):
    if os.path.splitext(str(path))[1].lower() != ".pdf":
        raise ValueError(f"RCTD figure output must use the .pdf extension: {path}")
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)


def load_rctd_weights(path):
    weights = pd.read_csv(path, index_col=0)
    weights.index = weights.index.astype(str)
    weights.columns = weights.columns.astype(str)

    if weights.empty or weights.shape[1] == 0:
        raise ValueError("RCTD weights must contain at least one spot and one cell type.")
    if not weights.index.is_unique or (weights.index.str.strip() == "").any():
        raise ValueError("RCTD weight row names must be unique, non-empty spot identifiers.")
    if not weights.columns.is_unique or (weights.columns.str.strip() == "").any():
        raise ValueError("RCTD weight columns must be unique, non-empty cell-type names.")

    try:
        weights = weights.apply(pd.to_numeric, errors="raise").astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("RCTD weights must be numeric.") from exc

    values = weights.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("RCTD weights contain non-finite values.")
    if np.any(values < -1e-8):
        raise ValueError("RCTD weights contain negative values.")

    weights = weights.clip(lower=0.0)
    row_sums = weights.sum(axis=1)
    if (~np.isfinite(row_sums)).any() or (row_sums <= 0).any():
        raise ValueError("RCTD weights contain rows with non-positive totals.")
    if not np.allclose(row_sums.to_numpy(), 1.0, rtol=1e-6, atol=1e-8):
        logger.warning("RCTD weights were not exactly row-normalized; normalizing them before use")
    return weights.div(row_sums, axis=0)


def load_doublet_results(path, expected_index):
    try:
        results = pd.read_csv(path, index_col=0)
    except pd.errors.EmptyDataError as exc:
        raise ValueError("RCTD doublet results file is empty.") from exc

    results.index = results.index.astype(str)
    if results.empty or not results.index.is_unique:
        raise ValueError("RCTD doublet results must contain unique spot identifiers.")

    required = {"spot_class", "first_type", "second_type"}
    missing = sorted(required.difference(results.columns))
    if missing:
        raise ValueError(f"RCTD doublet results are missing required columns: {missing}")
    if set(results.index) != set(expected_index):
        raise ValueError(
            "RCTD doublet results and weights must contain the same spot identifiers."
        )

    for column in ("spot_class", "first_type", "second_type"):
        results[column] = results[column].astype("string").str.strip()
        results.loc[results[column].eq(""), column] = pd.NA

    if results["spot_class"].isna().any():
        raise ValueError("RCTD doublet results contain missing spot_class values.")
    allowed_classes = {"singlet", "doublet_certain", "doublet_uncertain", "reject"}
    unexpected = sorted(set(results["spot_class"].astype(str)).difference(allowed_classes))
    if unexpected:
        raise ValueError(f"Unexpected RCTD spot_class value(s): {unexpected}")

    missing_first = results["spot_class"].ne("reject") & results["first_type"].isna()
    missing_second = results["spot_class"].eq("doublet_certain") & results["second_type"].isna()
    if missing_first.any():
        raise ValueError("Non-reject RCTD spots must have first_type predictions.")
    if missing_second.any():
        raise ValueError("doublet_certain RCTD spots must have second_type predictions.")
    return results


def get_spatial_coordinates(adata):
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"])
    elif all(column in adata.obs.columns for column in ("x", "y")):
        coords = adata.obs[["x", "y"]].to_numpy()
    else:
        raise ValueError("No spatial coordinates were found in obsm['spatial'] or obs[['x', 'y']].")

    if coords.ndim != 2 or coords.shape[0] != adata.n_obs or coords.shape[1] < 2:
        raise ValueError(
            f"Invalid spatial coordinate shape {coords.shape}; expected ({adata.n_obs}, >=2)."
        )
    coords = np.asarray(coords[:, :2], dtype=float)
    if not np.isfinite(coords).all():
        raise ValueError("Spatial coordinates contain non-finite values.")
    return coords


def align_weights_to_adata(adata, weights):
    adata.obs_names = adata.obs_names.astype(str)
    if not adata.obs_names.is_unique or (adata.obs_names.str.strip() == "").any():
        raise ValueError(
            "SpatialData table obs names must be unique, non-empty spot identifiers."
        )

    unknown_spots = weights.index.difference(adata.obs_names)
    if len(unknown_spots) > 0:
        raise ValueError(
            "RCTD weights contain spot identifiers that are absent from the source "
            f"SpatialData table; examples={list(unknown_spots[:5])}"
        )

    common = adata.obs_names.intersection(weights.index, sort=False)
    if len(common) == 0:
        raise ValueError(
            "No matching spot identifiers were found between SpatialData and RCTD weights. "
            f"Spatial examples={list(adata.obs_names[:3])}; RCTD examples={list(weights.index[:3])}"
        )

    aligned = pd.DataFrame(
        np.nan,
        index=pd.Index(adata.obs_names, name=adata.obs_names.name),
        columns=weights.columns,
        dtype=float,
    )
    aligned.loc[common, :] = weights.loc[common, :].to_numpy()
    adata.obsm["RCTD_weights"] = aligned

    omitted = adata.n_obs - len(common)
    logger.info("Matched %d/%d SpatialData spots to RCTD weights", len(common), adata.n_obs)
    if omitted:
        logger.warning(
            "%d input spot(s) were filtered before RCTD fitting; their RCTD weights remain NA",
            omitted,
        )
    return common


def plot_rctd_full_dotplot(
    adata,
    weights,
    cluster_col,
    sample_col,
    output_pdf,
    source_data_path,
    max_cell_types=30,
    enrichment_clip=2.5,
):
    """Use the cell2Location regional algorithm for normalized RCTD proportions."""

    if cluster_col not in adata.obs.columns:
        raise KeyError(
            f"Spatial cluster column '{cluster_col}' was not found in adata.obs; "
            "full-mode dotplot requires an existing unsupervised region label."
        )
    if sample_col is not None and sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column '{sample_col}' was not found in adata.obs.")
    _validate_pdf_path(output_pdf)
    os.makedirs(os.path.dirname(os.path.abspath(source_data_path)), exist_ok=True)

    metadata = adata.obs.loc[weights.index]
    cluster_text = metadata[cluster_col].astype("string").str.strip()
    valid = cluster_text.notna() & cluster_text.ne("")
    if sample_col is not None:
        sample_text = metadata[sample_col].astype("string").str.strip()
        valid &= sample_text.notna() & sample_text.ne("")
    if cluster_text.loc[valid].nunique() < 2:
        raise ValueError("At least two spatial regions are required for the RCTD full dotplot.")

    plot_adata = adata[weights.index].copy()
    abundance_key = "means_cell_abundance_w_sf"
    plot_adata.obsm[abundance_key] = weights.copy()
    model_metadata = dict(plot_adata.uns.get("mod", {}))
    model_metadata["factor_names"] = weights.columns.tolist()
    plot_adata.uns["mod"] = model_metadata

    result = plot_cell2location_dotplot(
        plot_adata,
        unsupervised_cluster_col=cluster_col,
        sample_col=sample_col,
        abundance_key=abundance_key,
        max_cell_types=max_cell_types,
        enrichment_clip=enrichment_clip,
        save_path=output_pdf,
        source_data_path=source_data_path,
        size_legend_title="Mean RCTD proportion",
        plot_title="Relative RCTD composition across spatial regions",
        x_axis_label="Unsupervised spatial region",
        method_label="RCTD",
    )

    source_data = result["source_data"].rename(
        columns={
            "mean_abundance": "mean_weight",
            "global_mean_abundance": "tissue_mean_weight",
        }
    )
    source_data = source_data.drop(columns=["abundance_key"], errors="ignore")
    source_data["weight_scale"] = "row-normalized RCTD proportion"
    source_data.to_csv(source_data_path, sep="\t", index=False)
    return source_data


def plot_rctd_doublet_proportion_dotplot(
    adata,
    cluster_col,
    sample_col,
    output_pdf,
    source_data_path,
    max_cell_types=30,
    enrichment_clip=2.5,
):
    """Summarize non-reject first-type proportions across spatial regions."""

    required = {"spot_class", "first_type", cluster_col}
    missing = sorted(required.difference(adata.obs.columns))
    if missing:
        raise KeyError(f"RCTD doublet proportion dotplot is missing obs columns: {missing}")
    if sample_col is not None and sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column '{sample_col}' was not found in adata.obs.")
    _validate_pdf_path(output_pdf)
    os.makedirs(os.path.dirname(os.path.abspath(source_data_path)), exist_ok=True)

    spot_class = adata.obs["spot_class"].astype("string").str.strip()
    first_type = adata.obs["first_type"].astype("string").str.strip()
    valid = spot_class.notna() & spot_class.ne("reject")
    valid &= first_type.notna() & first_type.ne("")
    if not valid.any():
        raise ValueError("No non-reject RCTD spots with first_type were available for plotting.")

    plot_adata = adata[valid.to_numpy()].copy()
    plotted_first_type = first_type.loc[valid].astype(str)
    factor_names = sorted(plotted_first_type.unique().tolist(), key=_natural_key)
    first_type_one_hot = pd.get_dummies(plotted_first_type).reindex(
        columns=factor_names,
        fill_value=0,
    ).astype(float)
    first_type_one_hot.index = plot_adata.obs_names.astype(str)

    abundance_key = "means_cell_abundance_w_sf"
    plot_adata.obsm[abundance_key] = first_type_one_hot
    model_metadata = dict(plot_adata.uns.get("mod", {}))
    model_metadata["factor_names"] = factor_names
    plot_adata.uns["mod"] = model_metadata

    result = plot_cell2location_dotplot(
        plot_adata,
        unsupervised_cluster_col=cluster_col,
        sample_col=sample_col,
        abundance_key=abundance_key,
        max_cell_types=max_cell_types,
        enrichment_clip=enrichment_clip,
        save_path=output_pdf,
        source_data_path=source_data_path,
        size_legend_title="Mean first-type proportion",
        plot_title="RCTD first-type composition across spatial regions",
        x_axis_label="Unsupervised spatial region",
        method_label="RCTD doublet",
    )

    source_data = result["source_data"].rename(
        columns={
            "mean_abundance": "mean_proportion",
            "global_mean_abundance": "tissue_mean_proportion",
        }
    )
    source_data = source_data.drop(columns=["abundance_key"], errors="ignore")
    source_data["denominator"] = "non-reject spots with first_type"
    source_data.to_csv(source_data_path, sep="\t", index=False)
    return source_data


def plot_rctd_spot_class_bar(adata, output_pdf, source_data_path):
    """Show counts and proportions for the official RCTD doublet classes."""

    if "spot_class" not in adata.obs.columns:
        raise KeyError("RCTD spot-class bar plot requires obs['spot_class'].")
    _validate_pdf_path(output_pdf)
    os.makedirs(os.path.dirname(os.path.abspath(source_data_path)), exist_ok=True)

    class_order = ["singlet", "doublet_certain", "doublet_uncertain", "reject"]
    display_labels = ["Singlet", "Certain\ndoublet", "Uncertain\ndoublet", "Reject"]
    spot_class = adata.obs["spot_class"].astype("string").str.strip().dropna()
    if spot_class.empty:
        raise ValueError("No fitted spots with RCTD spot_class were available for plotting.")

    counts = spot_class.value_counts().reindex(class_order, fill_value=0).astype(int)
    proportions = counts / counts.sum()
    source_data = pd.DataFrame(
        {
            "spot_class": class_order,
            "n_spots": counts.to_numpy(),
            "proportion": proportions.to_numpy(),
        }
    )
    source_data.to_csv(source_data_path, sep="\t", index=False)

    colors = ["#4C78A8", "#59A14F", "#F2CF5B", "#C96A64"]
    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig, ax = plt.subplots(figsize=(4.4, 3.0), constrained_layout=True)
        x_positions = np.arange(len(class_order))
        bars = ax.bar(
            x_positions,
            proportions.to_numpy() * 100,
            width=0.68,
            color=colors,
            edgecolor="none",
        )
        y_max = max(5.0, float(proportions.max() * 100) * 1.2)
        ax.set_ylim(0, min(105.0, y_max))
        for bar, count, proportion in zip(bars, counts, proportions):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + y_max * 0.025,
                f"{int(count):,}\n({proportion:.1%})",
                ha="center",
                va="bottom",
                fontsize=5.8,
            )
        ax.set_xticks(x_positions)
        ax.set_xticklabels(display_labels, fontsize=6)
        ax.set_ylabel("Spots (%)", fontsize=7)
        ax.set_title("RCTD spot classification", loc="left", fontsize=8, fontweight="semibold")
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.45)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", length=0, pad=3)
        ax.tick_params(axis="y", length=2, width=0.4, labelsize=5.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.savefig(
            output_pdf,
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.04,
            facecolor="white",
            metadata={"Title": "RCTD spot classification"},
        )
        plt.close(fig)

    logger.info("RCTD spot-class bar plot saved to %s", output_pdf)
    logger.info("RCTD spot-class source data saved to %s", source_data_path)
    return source_data


def _categorical_palette(labels):
    labels = sorted(set(labels), key=_natural_key)
    colors = []
    for cmap_name in ("tab20", "tab20b", "tab20c"):
        colors.extend(plt.get_cmap(cmap_name).colors)
    if len(labels) > len(colors):
        logger.warning(
            "More than %d RCTD cell types are displayed; additional categorical colors may be similar",
            len(colors),
        )
        colors.extend(
            plt.get_cmap("hsv")(
                np.linspace(0, 1, len(labels) - len(colors), endpoint=False)
            )
        )
    return {label: colors[index] for index, label in enumerate(labels)}


def _plot_assignment_panel(ax, coords, assignments, palette, title, point_size):
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        s=point_size,
        c="#D9D9D9",
        linewidths=0,
        alpha=0.65,
        rasterized=len(coords) > 20000,
    )
    valid = assignments.notna().to_numpy()
    if valid.any():
        assigned = assignments.loc[assignments.notna()].astype(str)
        ax.scatter(
            coords[valid, 0],
            coords[valid, 1],
            s=point_size * 1.25,
            c=assigned.map(palette).to_list(),
            linewidths=0,
            alpha=0.95,
            rasterized=valid.sum() > 20000,
        )
    else:
        ax.text(
            0.5,
            0.04,
            "No confident assignments",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=6,
            color="#666666",
        )
    ax.set_title(f"{title}\n(n = {int(valid.sum()):,})", fontsize=7.2, fontweight="semibold")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_rctd_doublet_spatial(adata, output_pdf, source_data_path):
    """Plot first type and only confident second type from RCTD doublet mode."""

    required = {"spot_class", "first_type", "second_type"}
    missing = sorted(required.difference(adata.obs.columns))
    if missing:
        raise KeyError(f"RCTD doublet visualization is missing obs columns: {missing}")
    _validate_pdf_path(output_pdf)
    os.makedirs(os.path.dirname(os.path.abspath(source_data_path)), exist_ok=True)

    coords = get_spatial_coordinates(adata)
    spot_class = adata.obs["spot_class"].astype("string").str.strip()
    first_raw = adata.obs["first_type"].astype("string").str.strip()
    second_raw = adata.obs["second_type"].astype("string").str.strip()
    first_assignment = first_raw.where(spot_class.notna() & spot_class.ne("reject"))
    second_assignment = second_raw.where(spot_class.eq("doublet_certain"))

    labels = set(first_assignment.dropna().astype(str)).union(
        second_assignment.dropna().astype(str)
    )
    palette = _categorical_palette(labels)
    point_size = float(np.clip(20000.0 / max(adata.n_obs, 1), 2.0, 10.0))

    source_data = pd.DataFrame(
        {
            "spot_id": adata.obs_names.astype(str),
            "x": coords[:, 0],
            "y": coords[:, 1],
            "spot_class": spot_class.to_numpy(),
            "first_type": first_raw.to_numpy(),
            "second_type": second_raw.to_numpy(),
            "first_type_plotted": first_assignment.to_numpy(),
            "second_type_plotted": second_assignment.to_numpy(),
        }
    )
    source_data.to_csv(source_data_path, sep="\t", index=False)

    n_labels = len(palette)
    legend_columns = min(8, max(1, int(np.ceil(n_labels / 6))))
    legend_rows = int(np.ceil(n_labels / legend_columns)) if n_labels else 0
    bottom_margin = min(0.45, 0.08 + 0.035 * legend_rows)
    fig_height = min(8.5, 3.35 + 0.18 * legend_rows)

    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, fig_height))
        _plot_assignment_panel(
            axes[0],
            coords,
            first_assignment,
            palette,
            "First predicted type (non-reject)",
            point_size,
        )
        _plot_assignment_panel(
            axes[1],
            coords,
            second_assignment,
            palette,
            "Second predicted type (certain doublets)",
            point_size,
        )
        fig.suptitle("RCTD doublet-mode spatial assignments", fontsize=8, fontweight="semibold")
        fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=bottom_margin, wspace=0.06)

        if palette:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=color,
                    markeredgecolor="none",
                    markersize=4,
                    label=label,
                )
                for label, color in palette.items()
            ]
            fig.legend(
                handles=handles,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.015),
                ncol=legend_columns,
                frameon=False,
                fontsize=5.3 if n_labels <= 30 else 4.8,
                columnspacing=0.9,
                handletextpad=0.35,
                title="Reference cell type",
                title_fontsize=5.8,
            )

        fig.savefig(
            output_pdf,
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.04,
            facecolor="white",
            metadata={"Title": "RCTD doublet-mode spatial assignments"},
        )
        plt.close(fig)

    logger.info("RCTD doublet-mode spatial assignments saved to %s", output_pdf)
    logger.info("RCTD doublet-mode spatial source data saved to %s", source_data_path)
    return first_assignment, second_assignment


def _store_results_in_obs(adata, results):
    results = results.reindex(adata.obs_names)
    for column in results.columns:
        values = results[column]
        if pd.api.types.is_numeric_dtype(values.dtype):
            adata.obs[column] = values.to_numpy()
        else:
            serializable = values.astype(object)
            serializable.loc[values.isna()] = np.nan
            adata.obs[column] = pd.Categorical(serializable)


def main():
    parser = argparse.ArgumentParser(description="Merge and visualize RCTD results in SpatialData")
    parser.add_argument(
        "--spatial_zarr",
        required=True,
        help="Source SpatialData zarr supplied internally by Snakemake",
    )
    parser.add_argument("--rctd_weights", required=True, help="Normalized RCTD weights CSV")
    parser.add_argument("--rctd_results", required=True, help="RCTD metadata/results CSV")
    parser.add_argument("--output_zarr", required=True, help="Output annotated SpatialData zarr")
    parser.add_argument("--mode", required=True, choices=("full", "doublet"), help="RCTD mode")
    parser.add_argument(
        "--run_type",
        required=True,
        help="Spatial technology type supplied by the workflow",
    )
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--output_plot", required=True, help="Output PDF path")
    parser.add_argument("--source_data", required=True, help="Figure source-data TSV path")
    parser.add_argument("--cluster_col", default="celltype", help="Unsupervised spatial region column")
    parser.add_argument("--sample_col", default="sample", help="Sample column for balanced summaries")
    parser.add_argument("--max_cell_types", type=int, default=30, help="Top cell types in regional dotplots; 0=all")
    parser.add_argument("--enrichment_clip", type=float, default=2.5, help="Symmetric color limit")
    args = parser.parse_args()

    sample_col = args.sample_col.strip() or None
    if args.max_cell_types < 0:
        parser.error("--max_cell_types must be >= 0; use 0 to display all cell types")
    if not np.isfinite(args.enrichment_clip) or args.enrichment_clip <= 0:
        parser.error("--enrichment_clip must be a positive finite number")

    log_step(logger, 1, 5, "loading SpatialData and RCTD weights")
    sdata = spd.read_zarr(args.spatial_zarr)
    table_keys = list(sdata.tables.keys())
    if not table_keys:
        raise ValueError("No tables were found in the SpatialData zarr input.")
    if len(table_keys) != 1:
        raise ValueError(
            "RCTD result write-back requires exactly one SpatialData table; "
            f"found: {table_keys}"
        )
    table_name = table_keys[0]
    adata = sdata[table_name]
    weights = load_rctd_weights(args.rctd_weights)
    common = align_weights_to_adata(adata, weights)
    weights_common = weights.loc[common].copy()
    adata.uns["RCTD"] = {
        "mode": args.mode,
        "weight_key": "RCTD_weights",
        "weight_scale": "row_normalized_proportion",
        "cell_types": weights.columns.astype(str).tolist(),
        "n_input_spots": int(adata.n_obs),
        "n_fitted_spots": int(len(common)),
    }

    log_step(logger, 2, 5, "merging mode-specific RCTD results")
    if args.mode == "doublet":
        _store_results_in_obs(
            adata,
            load_doublet_results(args.rctd_results, weights.index),
        )

    log_step(logger, 3, 5, "generating RCTD PDF visualization")
    if args.mode == "full":
        plot_rctd_full_dotplot(
            adata=adata,
            weights=weights_common,
            cluster_col=args.cluster_col,
            sample_col=sample_col,
            output_pdf=args.output_plot,
            source_data_path=args.source_data,
            max_cell_types=args.max_cell_types,
            enrichment_clip=args.enrichment_clip,
        )
    else:
        first_assignment, second_assignment = plot_rctd_doublet_spatial(
            adata=adata,
            output_pdf=args.output_plot,
            source_data_path=args.source_data,
        )
        figure_dir = os.path.dirname(os.path.abspath(args.output_plot))
        proportion_plot = os.path.join(
            figure_dir,
            f"{args.sample_id}_RCTD_doublet_proportion_dotplot.pdf",
        )
        proportion_source = os.path.join(
            figure_dir,
            f"{args.sample_id}_RCTD_doublet_proportion_dotplot_source.tsv",
        )
        spot_class_plot = os.path.join(
            figure_dir,
            f"{args.sample_id}_RCTD_spot_class_bar.pdf",
        )
        spot_class_source = os.path.join(
            figure_dir,
            f"{args.sample_id}_RCTD_spot_class_bar_source.tsv",
        )
        plot_rctd_doublet_proportion_dotplot(
            adata=adata,
            cluster_col=args.cluster_col,
            sample_col=sample_col,
            output_pdf=proportion_plot,
            source_data_path=proportion_source,
            max_cell_types=args.max_cell_types,
            enrichment_clip=args.enrichment_clip,
        )
        plot_rctd_spot_class_bar(
            adata=adata,
            output_pdf=spot_class_plot,
            source_data_path=spot_class_source,
        )
        adata.obs["RCTD_first_type"] = pd.Categorical(
            first_assignment.fillna("Unassigned").astype(str)
        )
        adata.obs["RCTD_second_type"] = pd.Categorical(
            second_assignment.fillna("Unassigned").astype(str)
        )

    log_step(logger, 4, 5, "saving RCTD results in SpatialData")
    sdata[table_name] = adata
    output_parent = os.path.dirname(os.path.abspath(args.output_zarr))
    os.makedirs(output_parent, exist_ok=True)

    if args.mode == "doublet":
        export_cluster_csv(
            sdata,
            args.run_type,
            output_parent,
            cell_id_col="cell_id",
            info_col="RCTD_first_type",
            sample_col=sample_col or "sample",
            sample_id=args.sample_id,
        )
    else:
        logger.info("Full mode stores continuous RCTD weights; no discrete cluster CSV is exported")

    sdata.write(args.output_zarr, overwrite=True)
    log_step(logger, 5, 5, f"RCTD outputs saved to {args.output_zarr}")
    logger.info("Figure saved to %s", args.output_plot)
    logger.info("RCTD visualization module completed")


if __name__ == "__main__":
    main()
