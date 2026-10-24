from __future__ import annotations

import pandas as pd
from spatialdata.transformations import get_transformation

from spatialsnake.workflow.function.stereoseq_spec import parse_stereoseq_input_spec

EMPTY_INPUT_SPECS = {None, "", False, "None", "False", "false", "NULL", "null"}
DEFAULT_REGION_KEY = "region"
DEFAULT_INSTANCE_KEY = "instance_id"


def resolve_region_image(sdata, table, region, system, image_type="hires"):
    """Match images to one annotated region, including shared coordinate systems."""
    images = [name for name in sdata.images
              if system in get_transformation(sdata.images[name], get_all=True)]
    if not images:
        return None
    attrs = table.uns.get("spatialdata_attrs", {})
    region_key = attrs.get("region_key")
    obs = table.obs
    if region_key in obs:
        obs = obs.loc[obs[region_key].astype(str).eq(str(region))]
    prefixes = [str(region)]
    for column in ("sample", "dataset_id"):
        if column in obs:
            prefixes.extend(obs[column].dropna().astype(str).unique())
    # Reader element names use these suffixes independently of sample IDs.
    for suffix in ("_polygons", "_cell_boundaries", "_cell_circles", "_cellbin_shapes"):
        if str(region).endswith(suffix):
            prefixes.append(str(region)[:-len(suffix)])
    matched = [name for name in images
               if any(name == prefix or name.startswith(prefix + "_") for prefix in prefixes)]
    if matched:
        images = matched
    elif "sample" in table.obs and table.obs["sample"].nunique() > 1:
        other_prefixes = set()
        for column in ("sample", "dataset_id"):
            if column in table.obs:
                other_prefixes.update(set(table.obs[column].dropna().astype(str)) - set(prefixes))
        if all(any(name == prefix or name.startswith(prefix + "_") for prefix in other_prefixes) for name in images):
            return None
        raise ValueError(f"Cannot match images to spatial region {region!r}")
    for keyword in (image_type, "hires", "morphology_focus", "morphology_mip", "HE_regist", "_z3"):
        preferred = [name for name in images if keyword and str(keyword) in name]
        if preferred:
            images = preferred
            break
    if len(images) != 1:
        raise ValueError(f"Multiple images match spatial region {region!r}: {images}; set image_type")
    return images[0]


def spatial_region_panels(sdata, table, image_type="hires", vis_mode="auto"):
    attrs = table.uns.get("spatialdata_attrs", {})
    region_key = attrs.get("region_key")
    regions = (pd.unique(table.obs[region_key].astype(str)).tolist()
               if region_key in table.obs else attrs.get("region", []))
    regions = [regions] if isinstance(regions, str) else list(regions)
    panels = []
    for region in regions:
        mode, element_name = resolve_visual_target(sdata, region, vis_mode)
        systems = sorted(get_transformation(sdata[element_name], get_all=True))
        if not systems:
            raise ValueError(f"Spatial region {region!r} has no coordinate system")
        system = next((name for name in systems if image_type and image_type in name), systems[0])
        image = resolve_region_image(sdata, table, region, system, image_type)
        panels.append((region, mode, element_name, system, image))
    return panels


def resolve_stereoseq_table_key(_input_spec, table_keys):
    table_keys = list(table_keys)
    if not table_keys:
        raise ValueError("No table found in SpatialData input.")

    non_analysis_tables = [key for key in table_keys if not str(key).startswith("analysis_")]
    fallback_tables = non_analysis_tables if non_analysis_tables else table_keys
    requested_tables = []

    if _input_spec not in EMPTY_INPUT_SPECS:
        parsed_specs = parse_stereoseq_input_spec(_input_spec)
        if parsed_specs:
            for item in parsed_specs:
                if item == "cellbin":
                    requested_tables.append("cellbin_table")
                elif item == "adjusted_cellbin":
                    requested_tables.append("adjusted_cellbin_table")
                else:
                    requested_tables.append(f"bin{int(item)}_table")

    for requested in requested_tables:
        if requested in fallback_tables:
            return requested
    if len(fallback_tables) == 1:
        return fallback_tables[0]
    raise ValueError(f"Select one Stereo-seq analysis table; requested {_input_spec!r}, available {fallback_tables}")


def extract_table_spatial_metadata(adata):
    spatial_attrs = adata.uns.get("spatialdata_attrs", {})
    region_name = spatial_attrs.get("region")
    if isinstance(region_name, (list, tuple)):
        region_name = region_name[0] if region_name else None
    return {
        "region_name": region_name,
        "region_key": spatial_attrs.get("region_key", DEFAULT_REGION_KEY),
        "instance_key": spatial_attrs.get("instance_key", DEFAULT_INSTANCE_KEY),
    }


def resolve_point_region_name(sdata, region_name):
    point_keys = set(getattr(sdata, "points", {}).keys())
    if not region_name:
        return None
    if region_name in point_keys:
        return region_name

    candidates = []
    if str(region_name).endswith("_shapes"):
        candidates.append(str(region_name)[: -len("_shapes")] + "_points")
    if str(region_name).endswith("_boundaries"):
        candidates.append(str(region_name)[: -len("_boundaries")] + "_points")
    if str(region_name).endswith("_circles"):
        candidates.append(str(region_name)[: -len("_circles")] + "_points")

    for candidate in candidates:
        if candidate in point_keys:
            return candidate
    return None


def resolve_visual_target(sdata, region_name, vis_mode=None):
    requested_mode = str(vis_mode or "auto").strip().lower()
    if requested_mode not in {"auto", "point", "shape"}:
        raise ValueError("vis_mode only supports 'auto', 'point', or 'shape'.")

    shape_keys = set(getattr(sdata, "shapes", {}).keys())
    has_shape = region_name in shape_keys if region_name else False
    point_region_name = resolve_point_region_name(sdata, region_name)
    has_point = point_region_name is not None

    if requested_mode == "auto":
        if has_shape:
            return "shape", region_name
        if has_point:
            return "point", point_region_name
        raise ValueError("No compatible shape or point element found for visualization.")

    if requested_mode == "shape":
        if not has_shape:
            raise ValueError(f"Requested shape visualization but shape region '{region_name}' was not found.")
        return "shape", region_name

    if not has_point:
        raise ValueError(f"Requested point visualization but no point region matched '{region_name}'.")
    return "point", point_region_name


def extract_point_coordinates_for_obs(sdata, adata, region_name=None, instance_key=None):
    point_region_name = resolve_point_region_name(sdata, region_name)
    if not point_region_name:
        return None

    points = sdata.points[point_region_name]
    if hasattr(points, "compute"):
        points = points.compute()
    if hasattr(points, "to_pandas") and not hasattr(points, "columns"):
        points = points.to_pandas()
    if not hasattr(points, "columns"):
        return None

    x_col = next((col for col in ("x", "center_x") if col in points.columns), None)
    y_col = next((col for col in ("y", "center_y") if col in points.columns), None)
    if x_col is None or y_col is None:
        return None

    coords = points[[x_col, y_col]].copy()
    coords.index = coords.index.astype(str)

    obs_names = adata.obs_names.astype(str)
    obs_ids = obs_names
    if instance_key and instance_key in adata.obs.columns:
        obs_ids = adata.obs[instance_key].astype(str)

    aligned = coords.reindex(pd.Index(obs_ids))
    aligned.index = obs_names
    if aligned[[x_col, y_col]].isna().all().all() and len(coords) == adata.n_obs:
        aligned = coords.iloc[: adata.n_obs].copy()
        aligned.index = obs_names

    aligned = aligned.dropna(subset=[x_col, y_col])
    if aligned.empty:
        return None
    aligned.columns = ["x", "y"]
    return aligned
