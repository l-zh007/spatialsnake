import os
import argparse
import numpy as np
import pandas as pd
import spatialdata as spd
import scanpy as sc
import matplotlib.pyplot as plt
from spatialdata import bounding_box_query
from spatialsnake.workflow.function.transform import zarr_to_h5ad
from spatialsnake.workflow.function.plot import plot_cell2location_dotplot
import spatialdata_plot as splt
import spatialdata_io as so

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", required=True)
parser.add_argument("--sample_id", required=True)
parser.add_argument("--output_zarr_path", required=True)
parser.add_argument("--type", required=True)
parser.add_argument("--image_type", required=True)
parser.add_argument("--image_slice", required=False)
parser.add_argument("--sample_cnt", type=int, required=True)
parser.add_argument("--coord", type=int, nargs=4, required=False)
parser.add_argument("--input_st", required=False)
parser.add_argument("--shape_type", required=False)
parser.add_argument("--celltype_col", required=False, default="celltype")
args = parser.parse_args()


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


def crop0(x, crs, x1, x2, y1, y2):
    return bounding_box_query(
        x,
        min_coordinate=[x1, y1],
        max_coordinate=[x2, y2],
        axes=("x", "y"),
        target_coordinate_system=crs
    )


def normalize_factor_names(factor_names):
    return [str(name).replace("+", "plus") for name in factor_names]


def ensure_factor_columns(adata, factor_names):
    normalized_names = normalize_factor_names(factor_names)
    adata.obs.columns = adata.obs.columns.astype(str).str.replace("+", "plus", regex=False)
    existing_cols = [col for col in normalized_names if col in adata.obs.columns]
    return normalized_names, existing_cols


def add_factor_columns_from_obsm(adata, factor_names):
    if not factor_names or "q05_cell_abundance_w_sf" not in adata.obsm:
        return []
    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if hasattr(abundance, "toarray"):
        abundance = abundance.toarray()
    if isinstance(abundance, pd.DataFrame):
        abundance_df = abundance.copy()
    else:
        abundance_df = pd.DataFrame(np.asarray(abundance), index=adata.obs_names)
    normalized_names = normalize_factor_names(factor_names)
    if abundance_df.shape[1] == len(normalized_names):
        abundance_df.columns = normalized_names
    else:
        abundance_df.columns = [f"cell_type_{i}" for i in range(abundance_df.shape[1])]
    abundance_df.index = adata.obs_names
    for column in abundance_df.columns:
        adata.obs[column] = abundance_df[column]
    return abundance_df.columns.tolist()


def prepare_scanpy_adata(zarr_path, table_adata, output_dir):
    converted = zarr_to_h5ad(
        [zarr_path],
        os.path.join(output_dir, "convert_data.h5ad"),
        save_image=True
    )
    if converted is None:
        return None
    converted.obs_names = converted.obs_names.astype(str)
    table_copy = table_adata.copy()
    table_copy.obs_names = table_copy.obs_names.astype(str)
    common_obs = converted.obs_names.intersection(table_copy.obs_names)
    if len(common_obs) == 0:
        return None
    converted = converted[common_obs].copy()
    table_copy = table_copy[common_obs].copy()
    for column in table_copy.obs.columns:
        converted.obs[column] = table_copy.obs[column]
    for key in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf", "q95_cell_abundance_w_sf", "stds_cell_abundance_w_sf"):
        if key in table_copy.obsm:
            converted.obsm[key] = table_copy.obsm[key]
    if "mod" in table_copy.uns:
        converted.uns["mod"] = table_copy.uns["mod"]
    return converted


def get_library_id(adata):
    if "library_id" in adata.obs and adata.obs["library_id"].notna().any():
        return str(adata.obs["library_id"].astype(str).iloc[0])
    spatial_uns = adata.uns.get("spatial", {})
    if isinstance(spatial_uns, dict) and spatial_uns:
        return next(iter(spatial_uns.keys()))
    return None


def get_img_key(adata, library_id):
    spatial_uns = adata.uns.get("spatial", {})
    if library_id is None or library_id not in spatial_uns:
        return None
    images = spatial_uns[library_id].get("images", {})
    if "hires" in images:
        return "hires"
    if images:
        return next(iter(images.keys()))
    return None


def plot_abundance_stats(adata, factor_names, save_dir="."):
    print("Generating cell abundance stats plot...")
    if "clusters" not in adata.obs and "q05_cell_abundance_w_sf" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="q05_cell_abundance_w_sf", n_neighbors=15)
        sc.tl.leiden(adata, key_added="clusters")
        sc.tl.umap(adata, min_dist=0.3, spread=1)
        sc.pl.umap(
            adata,
            color=["clusters"],
            size=30,
            color_map="RdPu",
            ncols=2,
            legend_loc="on data",
            legend_fontsize=20,
            show=False
        )
        plt.savefig(os.path.join(save_dir, "umap.png"), dpi=300, bbox_inches="tight")
        plt.close()
    valid_factors = [factor for factor in factor_names if factor in adata.obs.columns]
    if not valid_factors or "clusters" not in adata.obs.columns:
        print("No valid cell type columns found for plotting.")
        return
    df = adata.obs[["clusters"] + valid_factors].copy()
    cluster_means = df.groupby("clusters")[valid_factors].mean()
    cluster_percents = cluster_means.div(cluster_means.sum(axis=1), axis=0) * 100
    _, ax = plt.subplots(figsize=(12, 8))
    cluster_percents.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_title("Cell Type Abundance by Cluster")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Percentage Abundance (%)")
    ax.legend(title="Cell Types", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    save_path = os.path.join(save_dir, "cluster_abundance_stacked_bar.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    cluster_percents.to_csv(os.path.join(save_dir, "cluster_abundance_stats.csv"))
    print(f"Saved abundance plot to {save_path}")


def render_spatial_plots(sdata, shapes, images, systems, output_dir, image_slice, coord):
    x1, x2, y1, y2 = coord if coord else (None, None, None, None)
    total = min(len(shapes), len(images), len(systems))
    for i in range(total):
        sys = systems[0] if len(systems) < 2 else systems[i]
        _, axes = plt.subplots(2, 1, figsize=(20, 13))
        current_sdata = sdata
        if image_slice and coord:
            current_sdata = crop0(current_sdata, sys, x1, x2, y1, y2)
        current_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image", coordinate_systems=sys)
        current_sdata.pl.render_images(images[i]).pl.render_shapes(
            shapes[i], color="cellLoca_type"
        ).pl.show(ax=axes[1], coordinate_systems=sys)
        plt.savefig(f"{output_dir}/figure/{images[i]}_Most_rank_celltype.png", dpi=300, bbox_inches="tight")
        plt.close()


def plot_scanpy_spatial_panels(adata, factor_names, existing_cell_cols, output_dir):
    if adata is None or "spatial" not in adata.uns:
        print("Skip scanpy spatial panels because no usable AnnData spatial metadata was prepared.")
        return
    library_id = get_library_id(adata)
    img_key = get_img_key(adata, library_id)
    if img_key is None:
        print("Skip scanpy spatial panels because no image key was found in adata.uns['spatial'].")
        return
    if existing_cell_cols:
        sc.pl.spatial(
            adata,
            color=existing_cell_cols,
            show=False,
            img_key=img_key,
            library_id=library_id,
            size=1.3
        )
        plt.savefig(os.path.join(output_dir, "each_celltype.png"), dpi=300, bbox_inches="tight")
        plt.close()
    valid_factor_names = [name for name in factor_names if name in adata.obs.columns]
    if valid_factor_names:
        sc.pl.spatial(
            adata,
            color=valid_factor_names,
            show=False,
            img_key=img_key,
            library_id=library_id,
            size=1.3
        )
        plt.savefig(os.path.join(output_dir, "factor_namescelltype.png"), dpi=300, bbox_inches="tight")
        plt.close()


image_slice = parse_bool(args.image_slice)
output_dir = os.path.dirname(args.output_zarr_path)
figure_dir = os.path.join(output_dir, "figure")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)

concatenated_sdata = spd.read_zarr(args.input_dir)
table_key = next(iter(concatenated_sdata.tables.keys()))
adata_vis = concatenated_sdata[table_key]

factor_names = []
if "mod" in adata_vis.uns and "factor_names" in adata_vis.uns["mod"]:
    factor_names = list(adata_vis.uns["mod"]["factor_names"])

added_factor_cols = add_factor_columns_from_obsm(adata_vis, factor_names)
normalized_names, existing_cell_cols = ensure_factor_columns(adata_vis, factor_names)
if not existing_cell_cols and added_factor_cols:
    existing_cell_cols = [col for col in added_factor_cols if col in adata_vis.obs.columns]

if existing_cell_cols:
    abundance_matrix = adata_vis.obs[existing_cell_cols]
    adata_vis.obs["cellLoca_type"] = abundance_matrix.idxmax(axis=1)
else:
    adata_vis.obs["cellLoca_type"] = "unknown"

coord = tuple(args.coord) if args.coord else None
sample_cnt = args.sample_cnt
image_elements = list(concatenated_sdata.images.keys())
shape_elements = list(concatenated_sdata.shapes.keys())
valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
shape_count = len(shape_elements) // sample_cnt if sample_cnt else len(shape_elements)
shapes = []
images = []
systems = []

if len(valid_coord_systems) > 1:
    for image_name in image_elements:
        if args.image_type in image_name:
            images.append(image_name)
else:
    images = image_elements

if len(valid_coord_systems) > 1:
    for system_name in valid_coord_systems:
        if args.image_type in system_name:
            systems.append(system_name)
else:
    systems = valid_coord_systems

if shape_count > 1 and args.type != "visium" and args.shape_type:
    for shape_name in shape_elements:
        if args.shape_type in shape_name:
            shapes.append(shape_name)
else:
    shapes = shape_elements

if len(images) != len(shapes):
    raise ValueError(
        f"Check the spatial data to make sure that for every image there is a shape: "
        f"images={images}, shapes={shapes}"
    )

render_spatial_plots(
    concatenated_sdata,
    shapes,
    images,
    systems,
    output_dir,
    image_slice,
    coord
)

if (
    len(existing_cell_cols) > 0
    and args.celltype_col in adata_vis.obs.columns
    and "q05_cell_abundance_w_sf" in adata_vis.obsm
):
    plot_cell2location_dotplot(
        adata_vis,
        unsupervised_cluster_col=args.celltype_col,
        save_path=os.path.join(figure_dir, "cell2location_cluster_correlation_dotplot.png")
    )
else:
    print(
        f"Skip cell2location correlation dotplot because '{args.celltype_col}' was not found in obs, "
        f"no valid abundance columns were available, or 'q05_cell_abundance_w_sf' was missing in obsm."
    )

adata_for_colocation = prepare_scanpy_adata(args.input_dir, adata_vis, output_dir)
if adata_for_colocation is None:
    raise ValueError("Failed to prepare a Scanpy-compatible AnnData object for cell2location co-location.")

from cell2location import run_colocation

res_dict, adata_coloc = run_colocation(
    adata_for_colocation,
    model_name="CoLocatedGroupsSklearnNMF",
    train_args={
        "n_fact": np.arange(11, 13),
        "sample_name_col": "sample",
        "n_restarts": 16
    },
    model_kwargs={"alpha": 0.01, "init": "random", "nmf_kwd_args": {"tol": 0.000001}},
    export_args={
        "path": f"{output_dir}/CoLocatedComb/",
        "plot_histology": True,
        "scanpy_alpha_img": 1.0,
    }
)

concatenated_sdata[table_key] = adata_vis
concatenated_sdata.write(args.output_zarr_path)

if len(existing_cell_cols) > 0:
    plot_scanpy_spatial_panels(
        adata_coloc,
        normalize_factor_names(factor_names),
        existing_cell_cols,
        figure_dir
    )
    if "mod" in adata_vis.uns:
        plot_abundance_stats(adata_vis, normalize_factor_names(factor_names), figure_dir)

print(sda)
