import argparse
import json
import os
import spatialdata_plot
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import spatialdata as spd
from spatialdata.models import PointsModel
from spatialdata.transformations import Identity

from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.plot import cluster_proportion
from spatialsnake.workflow.function.stereoseq_selection import (
    extract_table_spatial_metadata,
    resolve_stereoseq_table_key,
    resolve_visual_target,
)
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

STEREOSEQ_TYPES = {"stereoseq", "StereoSeq", "Stereo-seq"}
logger = setup_logger("reannotation")


def parse_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    value_str = str(value).strip().lower()
    if value_str in ["true", "1", "yes", "y", "t"]:
        return True
    if value_str in ["false", "0", "no", "n", "f", "none", "null", ""]:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def select_elements(items, keyword, run_type):
    if run_type not in ["visium", "visium_HD", "visium_segment"]:
        return list(items)
    if not keyword:
        return list(items)
    selected = [item for item in items if keyword in item]
    return selected if selected else list(items)


def pick_table_key(sdata, run_type=None, input_spec=None):
    table_keys = list(getattr(sdata, "tables", {}).keys())
    if not table_keys:
        raise RuntimeError("No table found in SpatialData input.")
    if run_type in STEREOSEQ_TYPES:
        return resolve_stereoseq_table_key(input_spec, table_keys)
    return "table" if "table" in table_keys else table_keys[0]
def resolve_annotation_map(anno_data, sample_id):
    if sample_id in anno_data:
        return anno_data[sample_id]
    first_key = next(iter(anno_data.keys()))
    return anno_data[first_key]


def annotate_adata(adata, annotation_map):
    cluster_key = "recluster" if "recluster" in adata.obs.columns else "clusters"
    if cluster_key not in adata.obs.columns:
        raise ValueError("no recluster or clusters column found in table obs")
    mapped = adata.obs[cluster_key].astype(str).map(annotation_map)
    adata.obs["celltype"] = mapped.fillna("Unknown").astype("category")
    return adata


def ensure_obs_columns(adata, sample_id):
    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = adata.obs_names.astype(str)
    if "region" not in adata.obs.columns:
        adata.obs["region"] = str(sample_id)
    if "group" not in adata.obs.columns:
        adata.obs["group"] = adata.obs["region"].astype(str)
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["region"] = adata.obs["region"].astype(str).astype("category")
    adata.obs["group"] = adata.obs["group"].astype(str)
    return adata


def crop0(x, crs, x1, x2, y1, y2):
    return spd.bounding_box_query(
        x,
        min_coordinate=[x1, y1],
        max_coordinate=[x2, y2],
        axes=("x", "y"),
        target_coordinate_system=crs,
    )



def save_umap_plot(adata, cluster_key, output_png):
    sc.pl.umap(adata, color=[cluster_key], show=False)
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()


def render_spatial_plots(sdata, table_key, adata, run_type, vis_mode, point_size, image_slice, coord, output_dir, summary_png):
    table_meta = extract_table_spatial_metadata(adata)
    region_name = table_meta["region_name"]
    if isinstance(region_name, (list, tuple)):
        region_name = region_name[0] if region_name else None
    adata.obs_names = adata.obs_names.astype(str)
    sample_cnt = adata.obs["sample"].astype(str).nunique() if "sample" in adata.obs else 1
    image_elements = list(getattr(sdata, "images", {}).keys())
    shape_elements = list(getattr(sdata, "shapes", {}).keys())
    point_elements = list(getattr(sdata, "points", {}).keys())
    valid_coord_systems = sorted(getattr(sdata, "coordinate_systems", []))

    try:
        active_vis_mode, visual_region_name = resolve_visual_target(sdata, region_name, vis_mode)
    except Exception:
        active_vis_mode = "shape" if region_name in shape_elements else "point"
        visual_region_name = region_name

    point_visual_elements = []
    if active_vis_mode == "shape":
        visual_elements = [visual_region_name] if visual_region_name in shape_elements else []
    else:
        point_visual_elements = select_elements(point_elements, region_name, run_type)
        if not point_visual_elements:
            point_visual_elements = [visual_region_name] if visual_region_name in point_elements else point_elements
        visual_elements = list(point_visual_elements)
        keep_ids = adata.obs_names.astype(str)
        point_transformations = {coord_system: Identity() for coord_system in valid_coord_systems} if valid_coord_systems else {"global": Identity()}
        for current_visual_region in point_visual_elements:
            source_points = sdata.points[current_visual_region]
            points_df = source_points.compute() if hasattr(source_points, "compute") else source_points.copy()
            points_df = points_df.copy()
            points_df.index = points_df.index.astype(str)
            points_df = points_df.reindex(keep_ids).copy()
            points_df["celltype"] = (
                adata.obs["celltype"].astype(str).reindex(points_df.index).fillna("Unknown").astype("category")
            )
            points_df = points_df.dropna(subset=["x", "y"])
            sdata.points[current_visual_region] = PointsModel.parse(
                points_df.sort_index(),
                coordinates={"x": "x", "y": "y"},
                transformations=point_transformations,
            )

    images = select_elements(image_elements, "hires", run_type)
    systems = select_elements(valid_coord_systems, "hires", run_type)

    if len(images) == 0 or len(systems) == 0:
        if "spatial" in adata.obsm_keys():
            coords = adata.obsm["spatial"]
            codes = adata.obs["celltype"].astype("category").cat.codes
            plt.figure(figsize=(8, 8))
            plt.scatter(coords[:, 0], coords[:, 1], c=codes, s=5, cmap="tab20", linewidths=0)
            plt.gca().invert_yaxis()
            plt.savefig(summary_png, dpi=300, bbox_inches="tight")
            plt.close()
            return
        save_umap_plot(adata, "celltype", summary_png)
        return

    total = min(max(1, len(visual_elements) * sample_cnt), len(images), max(1, len(systems)))
    saved_summary = False
    for i in range(total):
        sys = systems[0] if len(systems) < 2 else systems[i]
        image_name = images[min(i, len(images) - 1)]
        title = image_name.replace("_hires_image", "")
        current_visual_region = visual_region_name
        if active_vis_mode == "point" and len(point_visual_elements) > 0:
            current_visual_region = point_visual_elements[0] if len(point_visual_elements) < 2 else point_visual_elements[min(i, len(point_visual_elements) - 1)]
        current_sdata = sdata
        if image_slice and coord is not None:
            current_sdata = crop0(current_sdata, sys, coord[0], coord[1], coord[2], coord[3])

        _, axes = plt.subplots(2, 1, figsize=(20, 13))
        axes = np.ravel(axes)
        current_sdata.pl.render_images(image_name).pl.show(
            ax=axes[0],
            title="image",
            coordinate_systems=sys,
        )
        if active_vis_mode == "shape":
            current_sdata.pl.render_shapes(
                current_visual_region,
                color="celltype",
                table_name=table_key,
                outline_alpha=0,
                fill_alpha=1,
            ).pl.show(
                ax=axes[1],
                coordinate_systems=sys,
                title=title,
            )
        else:
            current_sdata.pl.render_images(image_name).pl.render_points(
                current_visual_region,
                color="celltype",
                size=point_size,
                alpha=0.8,
            ).pl.show(
                ax=axes[1],
                coordinate_systems=sys,
                title=title,
            )
        per_image_png = os.path.join(output_dir, f"{image_name}_Clusters.png")
        plt.savefig(per_image_png, dpi=300, bbox_inches="tight")
        if not saved_summary:
            plt.savefig(summary_png, dpi=300, bbox_inches="tight")
            saved_summary = True
        plt.close()

    if not saved_summary and "spatial" in adata.obsm_keys():
        coords = adata.obsm["spatial"]
        codes = adata.obs["celltype"].astype("category").cat.codes
        plt.figure(figsize=(8, 8))
        plt.scatter(coords[:, 0], coords[:, 1], c=codes, s=5, cmap="tab20", linewidths=0)
        plt.gca().invert_yaxis()
        plt.savefig(summary_png, dpi=300, bbox_inches="tight")
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="reannotation")
    parser.add_argument("--input_dir", type=str, required=True)
    parser.add_argument("--sample_id", type=str, required=True)
    parser.add_argument("--output_zarr_path", type=str, required=True)
    parser.add_argument("--output_csv", type=str, required=True)
    parser.add_argument("--type", type=str, required=True)
    parser.add_argument("--anno_data", type=str, required=True)
    parser.add_argument("--image_slice", type=parse_bool, required=False)
    parser.add_argument("--vis_mode", type=str, required=False, default="auto")
    parser.add_argument("--point_size", type=float, required=False, default=2.5)
    parser.add_argument("--input_spec", type=str, required=False)
    parser.add_argument("--coord", type=int, nargs=4, required=False)
    args = parser.parse_args()

    anno_data = json.loads(args.anno_data)
    annotation_map = resolve_annotation_map(anno_data, args.sample_id)
    output_dir = os.path.dirname(args.output_csv)
    os.makedirs(output_dir, exist_ok=True)
    cluster_key = "celltype"
    umap_png = os.path.join(output_dir, "umap_recluster.png")
    spatial_png = os.path.join(output_dir, "spatial_clusters.png")

    log_step(logger, 1, 5, "loading annotation map and input data")
    sdata = spd.read_zarr(args.input_dir)
    table_key = pick_table_key(sdata, run_type=args.type, input_spec=args.input_spec)
    adata = sdata[table_key].copy()
    adata = ensure_obs_columns(adata, args.sample_id)
    adata = annotate_adata(adata, annotation_map)
    sdata[table_key] = adata
    logger.info(f"Reannotated {adata.n_obs} observations")
    log_step(logger, 2, 5, "exporting reannotation table")
    exported_files = export_cluster_csv(
        sdata,
        args.type,
        output_dir,
        cell_id_col="cell_id",
        info_col="celltype",
        sample_col="region",
        sample_id=args.sample_id
    )
    sdata.write(args.output_zarr_path, overwrite=True)

    log_step(logger, 3, 5, "saving cell type proportion plot")
    cluster_proportion(
        adata,
        sample_col="region",
        cluster_col="celltype",
        figsize=(12, 8),
        palette="tab20",
        sort_samples=None,
        sort_clusters=None,
        title="celltype_proportion",
        save_path=os.path.join(output_dir, "celltype_proportion.png"),
        dpi=300,
    )
    log_step(logger, 4, 5, "saving UMAP and spatial plots")
    save_umap_plot(adata, cluster_key, umap_png)
    render_spatial_plots(
        sdata=sdata,
        table_key=table_key,
        adata=adata,
        run_type=args.type,
        vis_mode=args.vis_mode,
        point_size=args.point_size,
        image_slice=args.image_slice,
        coord=args.coord,
        output_dir=output_dir,
        summary_png=spatial_png,
    )
    if len(exported_files) > 0:
        os.replace(exported_files[0], args.output_csv)
    log_step(logger, 5, 5, f"reannotation outputs saved to {output_dir}")
    logger.info("Reannotation module completed")


if __name__ == "__main__":
    main()
