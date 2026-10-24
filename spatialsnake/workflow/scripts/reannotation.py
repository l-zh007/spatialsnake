import argparse
import json
import os
import spatialdata_plot
import matplotlib.pyplot as plt
import scanpy as sc
import spatialdata as spd
from spatialdata.transformations import get_transformation

from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.plot import cluster_proportion
from spatialsnake.workflow.function.stereoseq_selection import (
    resolve_stereoseq_table_key,
    resolve_visual_target,
    resolve_region_image,
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


def sanitize_obs_index_name(adata):
    instance_key = adata.uns.get("spatialdata_attrs", {}).get("instance_key")
    if instance_key and instance_key in adata.obs and adata.obs_names.name == instance_key:
        safe_name = "__obs_index__"
        while safe_name in adata.obs:
            safe_name = f"_{safe_name}"
        adata.obs_names.name = safe_name
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
    adata.obs_names = adata.obs_names.astype(str)
    attrs = adata.uns.get("spatialdata_attrs", {})
    regions = attrs.get("region", [])
    regions = [regions] if isinstance(regions, str) else list(regions or [])
    images = select_elements(list(sdata.images), "hires", run_type)
    systems = sorted(sdata.coordinate_systems)
    panels = []
    for region in regions:
        try:
            mode, element_name = resolve_visual_target(sdata, region, vis_mode)
        except ValueError:
            continue
        element = sdata.shapes[element_name] if mode == "shape" else sdata.points[element_name]
        candidates = [system for system in systems if system in get_transformation(element, get_all=True)]
        if not candidates:
            continue
        system = next((system for system in candidates if "hires" in system.lower()), candidates[0])
        image = resolve_region_image(sdata, adata, region, system)
        panels.append((region, mode, element_name, system, image))

    if panels:
        rows = 2 if any(image for _, _, _, _, image in panels) else 1
        _, axes = plt.subplots(rows, len(panels), figsize=(10 * len(panels), 13 if rows == 2 else 8), squeeze=False)
        groups = adata.obs["celltype"].astype("category").cat.categories.tolist()
        instance_key = attrs.get("instance_key")
        region_key = attrs.get("region_key")
        for i, (region, mode, element_name, system, image) in enumerate(panels):
            current_sdata = crop0(sdata, system, *coord) if image_slice and coord is not None else sdata
            barcode_ax = axes[-1, i]
            if image:
                current_sdata.pl.render_images(image).pl.show(
                    ax=axes[0, i], title="image", coordinate_systems=system,
                )
            if mode == "shape":
                current_sdata.pl.render_shapes(
                    element_name, color="#9E9E9E", fill_alpha=1, outline_alpha=0,
                ).pl.render_shapes(
                    element_name, color="celltype", groups=groups, table_name=table_key,
                    na_color=None, fill_alpha=1, outline_alpha=0,
                ).pl.show(ax=barcode_ax, coordinate_systems=system, title=str(region))
            else:
                current_sdata.pl.render_points(
                    element_name, color="#9E9E9E", size=point_size, alpha=1,
                ).pl.show(ax=barcode_ax, coordinate_systems=system, title=str(region))
                region_obs = adata.obs
                if region_key in region_obs:
                    region_obs = region_obs[region_obs[region_key].astype(str) == str(region)]
                label_ids = region_obs[instance_key].astype(str) if instance_key in region_obs else region_obs.index.astype(str)
                labels_by_id = region_obs["celltype"].copy()
                labels_by_id.index = label_ids
                points = spd.transform(current_sdata.points[element_name], to_coordinate_system=system).compute()
                points.index = points.index.astype(str)
                labels = labels_by_id.reindex(points.index).astype("category").cat.remove_unused_categories()
                colors = plt.get_cmap("tab20", max(1, len(labels.cat.categories)))
                for j, label in enumerate(labels.cat.categories):
                    selected = labels == label
                    barcode_ax.scatter(
                        points.loc[selected, "x"], points.loc[selected, "y"],
                        s=point_size, color=colors(j), label=label,
                    )
                if len(labels.cat.categories):
                    barcode_ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
        plt.tight_layout()
        plt.savefig(summary_png, dpi=300, bbox_inches="tight")
        plt.close()
        return

    if "spatial" in adata.obsm_keys():
        coords = adata.obsm["spatial"]
        codes = adata.obs["celltype"].astype("category").cat.codes
        plt.figure(figsize=(8, 8))
        plt.scatter(coords[:, 0], coords[:, 1], c=codes, s=5, cmap="tab20", linewidths=0)
        plt.gca().invert_yaxis()
        plt.savefig(summary_png, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        save_umap_plot(adata, "celltype", summary_png)

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
    umap_png = os.path.join(output_dir, "umap_annotations.png")
    spatial_png = os.path.join(output_dir, "spatial_annotations.png")

    log_step(logger, 1, 5, "loading annotation map and input data")
    sdata = spd.read_zarr(args.input_dir)
    table_key = pick_table_key(sdata, run_type=args.type, input_spec=args.input_spec)
    adata = sdata[table_key].copy()
    adata = sanitize_obs_index_name(ensure_obs_columns(adata, args.sample_id))
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
