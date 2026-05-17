import argparse
import os
import pandas as pd
import spatialdata as spd
import matplotlib.pyplot as plt
import numpy as np
import warnings
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv

warnings.filterwarnings("ignore")


def get_spatial_coordinates(adata):
    if "spatial" in adata.obsm:
        return np.asarray(adata.obsm["spatial"])
    if all(c in adata.obs.columns for c in ["x", "y"]):
        return adata.obs[["x", "y"]].to_numpy()
    return None


def plot_full_mode(adata, sample_id, output_plot):
    coords = get_spatial_coordinates(adata)
    if coords is None:
        raise ValueError("No valid spatial coordinates found for full-mode visualization.")

    labels = adata.obs["RCTD_max_celltype"].fillna("Unknown").astype(str)
    label_cat = pd.Categorical(labels)
    max_weight = adata.obs["RCTD_max_weight"].fillna(0).astype(float).to_numpy()

    n_labels = max(len(label_cat.categories), 1)
    cmap = plt.cm.get_cmap("tab20", n_labels)
    label_colors = label_cat.codes

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    axes[0].scatter(
        coords[:, 0],
        coords[:, 1],
        c=label_colors,
        s=10 + 25 * max_weight,
        cmap=cmap,
        linewidths=0,
        alpha=0.95
    )
    axes[0].set_title(f"{sample_id} - RCTD dominant cell type", fontsize=14, fontweight="bold")
    axes[0].set_aspect("equal")
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    legend_handles = []
    for idx, label in enumerate(label_cat.categories):
        legend_handles.append(
            plt.Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                color=cmap(idx),
                label=label,
                markersize=6
            )
        )
    axes[0].legend(
        handles=legend_handles,
        title="Dominant cell type",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False
    )

    scatter_weight = axes[1].scatter(
        coords[:, 0],
        coords[:, 1],
        c=max_weight,
        s=14,
        cmap="magma",
        linewidths=0,
        alpha=0.95
    )
    axes[1].set_title(f"{sample_id} - RCTD dominant proportion", fontsize=14, fontweight="bold")
    axes[1].set_aspect("equal")
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    cbar = fig.colorbar(scatter_weight, ax=axes[1], fraction=0.046, pad=0.04)
    cbar.set_label("Max normalized proportion", fontsize=11)

    for ax in axes:
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.tight_layout()
    fig.savefig(output_plot, bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_non_full_mode(sdata, adata, sample_id, output_plot):
    image_elements = list(sdata.images.keys())
    shape_elements = list(sdata.shapes.keys())
    valid_coord_systems = sorted(sdata.coordinate_systems)

    shapes_to_plot = []
    images_to_plot = []
    systems_to_plot = []

    hires_images = [img for img in image_elements if "hires" in img]
    if hires_images:
        images_to_plot = hires_images
    elif image_elements:
        images_to_plot = [image_elements[0]]

    if shape_elements:
        shapes_to_plot = [shape_elements[0]]

    if valid_coord_systems:
        systems_to_plot = [valid_coord_systems[0]]

    img_key = images_to_plot[0] if images_to_plot else None
    shape_key = shapes_to_plot[0] if shapes_to_plot else None
    sys_key = systems_to_plot[0] if systems_to_plot else None

    viz_col = "first_type" if "first_type" in adata.obs.columns else "RCTD_max_celltype"
    title = f"{sample_id} - RCTD Top Cell Type"

    fig, ax = plt.subplots(figsize=(10, 10))
    plotted = False

    if img_key:
        try:
            if sys_key:
                sdata.pl.render_images(img_key).pl.show(ax=ax, coordinate_systems=sys_key, title="")
            else:
                sdata.pl.render_images(img_key).pl.show(ax=ax, title="")
            plotted = True
        except Exception as e:
            print(f"Warning: Failed to render image {img_key}: {e}")

    if shape_key:
        try:
            if sys_key:
                sdata.pl.render_shapes(shape_key, color=viz_col).pl.show(
                    ax=ax,
                    coordinate_systems=sys_key,
                    title=title
                )
            else:
                sdata.pl.render_shapes(shape_key, color=viz_col).pl.show(ax=ax, title=title)
            plotted = True
        except Exception as e:
            print(f"Warning: Failed to render shapes {shape_key}: {e}")

    if not plotted:
        coords = get_spatial_coordinates(adata)
        if coords is None:
            raise ValueError("No valid coordinate system or coordinates found for plotting.")
        labels = pd.Categorical(adata.obs[viz_col].fillna("Unknown"))
        ax.scatter(coords[:, 0], coords[:, 1], c=labels.codes, s=5, cmap="tab20")
        ax.set_title(title)
        ax.set_aspect("equal")

    fig.savefig(output_plot, bbox_inches="tight", dpi=300)
    plt.close(fig)



def main():
    parser = argparse.ArgumentParser(description="Visualize RCTD results on SpatialData zarr")
    parser.add_argument("--zarr_input", type=str, required=True, help="Path to input SpatialData zarr")
    parser.add_argument("--rctd_weights", type=str, required=True, help="Path to RCTD weights CSV")
    parser.add_argument("--rctd_results", type=str, required=True, help="Path to RCTD metadata/results CSV")
    parser.add_argument("--output_zarr", type=str, required=True, help="Path to output annotated SpatialData zarr")
    parser.add_argument("--mode", type=str, required=True, help="RCTD mode")
    parser.add_argument("--sample_id", type=str, required=True, help="Sample ID")
    parser.add_argument("--output_plot", type=str, required=True, help="Path to save the visualization plot")
    args = parser.parse_args()

    # 1. Load SpatialData
    print(f"Loading SpatialData from: {args.zarr_input}")
    sdata = spd.read_zarr(args.zarr_input)
    
    # Assuming standard table structure
    # Try to find the table (AnnData)
    table_keys = list(sdata.tables.keys())
    if not table_keys:
        raise ValueError("No tables found in SpatialData zarr")
    
    table_name = table_keys[0] # Default to the first table
    adata = sdata[table_name]
    print(f"Using table: {table_name}")

    # 2. Load RCTD Results
    print(f"Loading RCTD weights from: {args.rctd_weights}")
    rctd_weights_df = pd.read_csv(args.rctd_weights, index_col=0)
    rctd_weights_df.index = rctd_weights_df.index.astype(str)
    print(f"Loading RCTD metadata/results from: {args.rctd_results}")
    rctd_results_df = pd.read_csv(args.rctd_results, index_col=0)
    rctd_results_df.index = rctd_results_df.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    
    # 3. Merge Results into AnnData
    common_indices = adata.obs_names.intersection(rctd_weights_df.index)
    if len(common_indices) == 0:
        print("Warning: No common indices found between SpatialData and RCTD results!")
        print("AnnData indices example:", adata.obs_names[:5])
        print("RCTD indices example:", rctd_weights_df.index[:5])
    else:
        print(f"Found {len(common_indices)} common spots/cells.")

    if rctd_weights_df.shape[1] == 0:
        raise ValueError("RCTD weights contain no cell type columns.")

    rctd_weights_reindexed = rctd_weights_df.reindex(adata.obs_names)
    rctd_results_reindexed = rctd_results_df.reindex(adata.obs_names)
    
    rctd_weights = rctd_weights_reindexed.fillna(0)
    adata.obsm["RCTD_weights"] = rctd_weights

    adata.obs["RCTD_max_celltype"] = rctd_weights.idxmax(axis=1).astype(str)
    adata.obs["RCTD_max_weight"] = rctd_weights.max(axis=1)

    for column in rctd_results_reindexed.columns:
        adata.obs[column] = rctd_results_reindexed[column]

    print(f"Saving annotated SpatialData to: {args.output_zarr}")
    os.makedirs(os.path.dirname(args.output_zarr), exist_ok=True)

    mode = str(args.mode).lower()
    info_col = "RCTD_max_celltype" if mode == "full" else "first_type"
    if info_col not in adata.obs.columns:
        info_col = "RCTD_max_celltype"

    export_dir = os.path.dirname(args.output_zarr) or "."
    data_type = "xenium" if "xenium" in f"{args.zarr_input}{args.output_zarr}".lower() else "visium"
    export_cluster_csv(
        sdata,
        data_type,
        export_dir,
        cell_id_col="cell_id",
        info_col=info_col,
        sample_col="sample",
        sample_id=args.sample_id
    )

    sdata.write(args.output_zarr, overwrite=True)

    print(f"Generating plots to: {args.output_plot}")

    if mode == "full":
        plot_full_mode(adata, args.sample_id, args.output_plot)
    else:
        plot_non_full_mode(sdata, adata, args.sample_id, args.output_plot)

    print(f"Plot saved to {args.output_plot}")

if __name__ == "__main__":
    main()
