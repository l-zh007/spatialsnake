from __future__ import annotations

import pandas as pd

from spatialsnake.workflow.function.stereoseq_spec import parse_stereoseq_input_spec

EMPTY_INPUT_SPECS = {None, "", False, "None", "False", "false", "NULL", "null"}
DEFAULT_REGION_KEY = "region"
DEFAULT_INSTANCE_KEY = "instance_id"


def resolve_stereoseq_table_key(_input_spec, table_keys):
    table_keys = list(table_keys)
    if not table_keys:
        raise ValueError("No table found in SpatialData input.")

    non_analysis_tables = [key for key in table_keys if not str(key).startswith("analysis_")]
    fallback_tables = non_analysis_tables if non_analysis_tables else table_keys
    requested_tables = []

    if _input_spec not in EMPTY_INPUT_SPECS:
        try:
            parsed_specs = parse_stereoseq_input_spec(_input_spec)
        except ValueError:
            parsed_specs = None
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
    for requested in requested_tables:
        requested_prefix = requested.replace("_table", "")
        for table_key in fallback_tables:
            if requested_prefix in str(table_key):
                return table_key
    if "table" in fallback_tables:
        return "table"
    return fallback_tables[0]


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
