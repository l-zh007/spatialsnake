"""Split AnnData/SpatialData objects without breaking downstream SpatialData use.

The default behaviour deliberately subsets only the selected table.  This is
important for the historical ``useful_tool --option=splitting`` workflow:
reclustering and reannotation consume the smaller table while still expecting
the original images, labels, shapes, points and transformations to be present.
Spatial elements are filtered only when ``--subset_mode spatial`` is requested.
"""

import argparse
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import spatialdata as spd
from shapely.geometry import Point, Polygon

from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("useful_splitting")


def build_parser():
    parser = argparse.ArgumentParser(description="Split an AnnData or SpatialData object")
    parser.add_argument("--INPUT_FILE", type=str, required=True, help="Input .h5ad file or SpatialData Zarr directory")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument("--split_by", type=str, required=True, help="obs column, ROI, or image")
    parser.add_argument("--barcodes", type=str, default="", help="Comma-combined or pipe-separated labels")
    # These arguments remain required for compatibility with the existing ToolRunner.
    parser.add_argument("--max_x", type=float, required=True, help="Image crop maximum x")
    parser.add_argument("--min_x", type=float, required=True, help="Image crop minimum x")
    parser.add_argument("--max_y", type=float, required=True, help="Image crop maximum y")
    parser.add_argument("--min_y", type=float, required=True, help="Image crop minimum y")
    parser.add_argument("--shape_elements", type=str, required=True, help="Legacy coordinate-system option")
    parser.add_argument("--roi_csv", type=str, default="", help="ROI CSV file or directory")
    parser.add_argument("--table_key", type=str, default="", help="Table to split when a SpatialData object has multiple tables")
    parser.add_argument("--subset_mode", choices=("table", "spatial"), default="table",
                        help="Subset only the table (default) or also associated spatial elements")
    parser.add_argument("--coordinate_system", type=str, default="", help="Coordinate system for image/polygon ROI coordinates")
    parser.add_argument("--roi_label_col", type=str, default="", help="Column containing ROI/category labels")
    parser.add_argument("--roi_region", type=str, default="", help="Limit ROI IDs to one integrated region")
    parser.add_argument("--annotation_format", choices=("auto", "loupe", "xenium", "both", "none"),
                        default="auto", help="Companion annotation CSV format")
    parser.add_argument("--annotation_cols", type=str, default="", help="Additional comma-separated obs columns to export")
    return parser


def safe_name(value):
    value = re.sub(r"[^\w.\-]+", "_", str(value).strip())
    return value.strip("_") or "unnamed"


def parse_barcode_groups(barcodes):
    """Parse commas as one combined subset and pipes as separate subsets."""
    text = "" if barcodes is None else str(barcodes).strip()
    if not text or text.lower() in {"none", "null", "na"}:
        return []
    groups = []
    for group_text in text.split("|"):
        values = [value.strip() for value in group_text.split(",") if value.strip()]
        if values:
            groups.append(values)
    return groups


def _normalise_column(value):
    return re.sub(r"[\s_\-]+", "", str(value).replace("\ufeff", "").strip().lower())


def _optional_text(value):
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in {"", "none", "null", "na"} else text


def _spatialdata_attrs(adata):
    attrs = adata.uns.get("spatialdata_attrs", {})
    return attrs if isinstance(attrs, dict) else {}


def _table_keys(adata):
    attrs = _spatialdata_attrs(adata)
    return attrs.get("region_key"), attrs.get("instance_key")


def prepare_subset_table(adata, mask):
    """Copy an AnnData subset and update only the declared region list."""
    subset = adata[mask, :].copy()
    if subset.n_obs == 0:
        raise ValueError("the requested split contains no observations")
    attrs = dict(_spatialdata_attrs(subset))
    region_key = attrs.get("region_key")
    if region_key and region_key in subset.obs.columns:
        attrs["region"] = pd.unique(subset.obs[region_key]).tolist()
        subset.uns["spatialdata_attrs"] = attrs
    return subset


def _copy_sdata_with_table(sdata, table_key, table):
    """Create a new container while retaining all non-selected elements."""
    tables = dict(sdata.tables)
    tables[table_key] = table
    return spd.SpatialData(
        images=dict(sdata.images),
        labels=dict(sdata.labels),
        points=dict(sdata.points),
        shapes=dict(sdata.shapes),
        tables=tables,
        attrs=dict(sdata.attrs),
    )


def _merge_rasters(filtered, source):
    """Keep raster context when instance-linked vectors are spatially filtered."""
    return spd.SpatialData(
        images=dict(source.images),
        labels=dict(source.labels),
        points=dict(filtered.points),
        shapes=dict(filtered.shapes),
        tables=dict(filtered.tables),
        attrs=dict(source.attrs),
    )


def _sdata_for_table_subset(sdata, table_key, table, subset_mode):
    if subset_mode == "table":
        return _copy_sdata_with_table(sdata, table_key, table)
    # match_sdata_to_table uses the table's original region_key/instance_key and
    # is consequently safer than guessing an element index or a cell_id column.
    filtered = spd.match_sdata_to_table(sdata, table_name=table_key, table=table, how="right")
    return _merge_rasters(filtered, sdata)


def _detect_input(path):
    path_obj = Path(path)
    if path_obj.is_file() and path_obj.suffix.lower() == ".h5ad":
        return "h5ad"
    if path_obj.is_dir():
        # SpatialData stores may use either Zarr v2 or v3 metadata and do not
        # have to carry a .zarr suffix.
        markers = (".zgroup", ".zattrs", "zarr.json")
        if path_obj.suffix.lower() == ".zarr" or any((path_obj / marker).exists() for marker in markers):
            return "zarr"
    raise ValueError(f"unsupported input {path!r}; expected .h5ad or a SpatialData Zarr directory")


def _select_table(sdata, requested):
    keys = list(sdata.tables.keys())
    if not keys:
        raise ValueError("SpatialData input contains no tables")
    if requested:
        if requested not in keys:
            raise ValueError(f"table_key {requested!r} not found; available tables: {', '.join(keys)}")
        return requested
    if len(keys) == 1:
        return keys[0]
    if "table" in keys:
        logger.info("Multiple tables found; selecting conventional table key 'table'")
        return "table"
    raise ValueError(f"multiple tables found ({', '.join(keys)}); specify --table_key")


def resolve_coordinate_system(sdata, requested="", legacy="", roi_region=""):
    systems = sorted(sdata.coordinate_systems)
    for raw_candidate in (requested, legacy, roi_region):
        candidate = _optional_text(raw_candidate)
        if candidate:
            if candidate not in systems:
                raise ValueError(f"coordinate system {candidate!r} not found; available: {', '.join(systems)}")
            return candidate
    if len(systems) == 1:
        return systems[0]
    raise ValueError(f"multiple coordinate systems found ({', '.join(systems)}); specify --coordinate_system")


def resolve_systems(sdata):
    return sorted(sdata.coordinate_systems)


def _platform_hint(sdata, adata, roi_hint=None):
    if roi_hint in {"loupe", "xenium"}:
        return roi_hint
    names = " ".join(
        list(getattr(sdata, "images", {}).keys())
        + list(getattr(sdata, "shapes", {}).keys())
        + list(getattr(sdata, "points", {}).keys())
    ).lower() if sdata is not None else ""
    obs_norm = {_normalise_column(column) for column in adata.obs.columns}
    if {"arrayrow", "arraycol"}.issubset(obs_norm) or "intissue" in obs_norm:
        return "loupe"
    if "barcode" in obs_norm or any(token in names for token in ("visium", "spot")):
        return "loupe"
    if (
        {"transcriptcounts", "cellarea"}.issubset(obs_norm)
        or "segmentationmethod" in obs_norm
        or any(token in names for token in ("xenium", "transcript", "cell_boundar"))
    ):
        return "xenium"
    return None


def _annotation_ids(adata, platform):
    region_key, instance_key = _table_keys(adata)
    norm_cols = {_normalise_column(column): column for column in adata.obs.columns}
    if platform == "loupe" and "barcode" in norm_cols:
        return adata.obs[norm_cols["barcode"]].astype(str), region_key
    if instance_key and instance_key in adata.obs.columns:
        return adata.obs[instance_key].astype(str), region_key
    if "cellid" in norm_cols:
        return adata.obs[norm_cols["cellid"]].astype(str), region_key
    return pd.Series(adata.obs_names.astype(str), index=adata.obs_names), region_key


def export_annotation_csvs(adata, output_path, split_by, annotation_format, annotation_cols,
                           sdata=None, roi_name=None, roi_hint=None):
    if annotation_format == "none":
        return
    if annotation_format == "auto":
        detected = _platform_hint(sdata, adata, roi_hint)
        if detected is None:
            logger.warning("Could not infer Loupe or Xenium annotation format; use --annotation_format to export CSV")
            return
        formats = [detected]
    elif annotation_format == "both":
        formats = ["loupe", "xenium"]
    else:
        formats = [annotation_format]

    requested_cols = []
    if split_by in adata.obs.columns:
        requested_cols.append(split_by)
    for column in str(annotation_cols or "").split(","):
        column = column.strip()
        if column and column not in requested_cols:
            requested_cols.append(column)
    missing_cols = [column for column in requested_cols if column not in adata.obs.columns]
    for column in missing_cols:
        logger.warning("Annotation column %r not found and will not be exported", column)
    requested_cols = [column for column in requested_cols if column in adata.obs.columns]
    if roi_name is not None:
        values = {"ROI": pd.Series(str(roi_name), index=adata.obs_names)}
    else:
        values = {}
    values.update({column: adata.obs[column].astype(str) for column in requested_cols})
    if not values:
        logger.warning("No annotation columns are available for CSV export")
        return

    output_stem = os.path.splitext(output_path)[0]
    for platform in formats:
        ids, region_key = _annotation_ids(adata, platform)
        if region_key and region_key in adata.obs.columns:
            regions = pd.unique(adata.obs[region_key].astype(str)).tolist()
        else:
            regions = [None]
        multiple_regions = len(regions) > 1
        for region in regions:
            mask = pd.Series(True, index=adata.obs_names)
            if region is not None:
                mask = adata.obs[region_key].astype(str).eq(region)
            suffix = f"_{safe_name(region)}" if multiple_regions else ""
            if platform == "loupe":
                frame = pd.DataFrame({"barcode": ids.loc[mask].to_numpy()}, index=adata.obs_names[mask])
                for column, series in values.items():
                    frame[column] = series.loc[mask].to_numpy()
                frame.to_csv(f"{output_stem}{suffix}_loupe.csv", index=False)
            else:
                for column, series in values.items():
                    frame = pd.DataFrame({
                        "cell_id": ids.loc[mask].to_numpy(),
                        "group": series.loc[mask].to_numpy(),
                    })
                    frame.to_csv(f"{output_stem}{suffix}_xenium_{safe_name(column)}.csv", index=False)


def _write_subset(subset, output_path, input_kind, sdata, table_key, args, split_by,
                  roi_name=None, roi_hint=None, polygon=None, coordinate_system=None):
    if input_kind == "h5ad":
        subset.write_h5ad(output_path)
        output_sdata = None
    else:
        if polygon is not None and args.subset_mode == "spatial":
            queried = spd.polygon_query(
                sdata,
                polygon=polygon,
                target_coordinate_system=coordinate_system,
                filter_table=True,
                clip=False,
            )
            if queried is None:
                raise ValueError(f"polygon ROI {roi_name!r} did not intersect the SpatialData object")
            filtered = spd.match_sdata_to_table(queried, table_name=table_key, table=subset, how="right")
            output_sdata = _merge_rasters(filtered, queried)
        else:
            output_sdata = _sdata_for_table_subset(sdata, table_key, subset, args.subset_mode)
        output_sdata.write(output_path, overwrite=True)
    export_annotation_csvs(
        subset,
        output_path,
        split_by,
        args.annotation_format,
        args.annotation_cols,
        sdata=output_sdata if input_kind == "zarr" else None,
        roi_name=roi_name,
        roi_hint=roi_hint,
    )


def split_metadata(adata, input_kind, sdata, table_key, args):
    split_by = args.split_by
    if split_by not in adata.obs.columns:
        raise ValueError(f"split_by column {split_by!r} not found in table obs")
    obs_values = adata.obs[split_by].astype(str)
    groups = parse_barcode_groups(args.barcodes)
    if not groups:
        groups = [[str(value)] for value in pd.unique(obs_values)]
        individual = True
    else:
        individual = False
    available = set(pd.unique(obs_values))
    for values in groups:
        missing = [value for value in values if value not in available]
        if missing:
            raise ValueError(f"values not found in {split_by}: {', '.join(missing)}")
        subset = prepare_subset_table(adata, obs_values.isin(values).to_numpy())
        if individual:
            filename = f"cluster_{safe_name(values[0])}.{input_kind}"
        else:
            filename = f"{safe_name(split_by)}_selected_{safe_name('_'.join(values))}.{input_kind}"
        output_path = os.path.join(args.output_dir, filename)
        _write_subset(subset, output_path, input_kind, sdata, table_key, args, split_by)
        logger.info("Wrote %s observations to %s", subset.n_obs, output_path)


def split_sample_or_group(adata, input_kind, sdata, table_key, args):
    if input_kind == "h5ad":
        split_metadata(adata, input_kind, sdata, table_key, args)
        return
    split_by = args.split_by
    if split_by in {"sample", "samples", "region"}:
        obs_values = adata.obs[split_by].astype(str)
        groups = parse_barcode_groups(args.barcodes)
        individual = not groups
        if individual:
            groups = [[str(value)] for value in pd.unique(obs_values)]
        available = set(pd.unique(obs_values))
        systems = set(sdata.coordinate_systems)
        for values in groups:
            missing = [value for value in values if value not in available]
            if missing:
                raise ValueError(f"values not found in {split_by}: {', '.join(missing)}")
            subset_table = prepare_subset_table(adata, obs_values.isin(values).to_numpy())
            candidate_systems = [value for value in values if value in systems]
            for source_col in ("sample", _table_keys(adata)[0]):
                if source_col and source_col in subset_table.obs.columns:
                    for value in pd.unique(subset_table.obs[source_col].astype(str)):
                        if value in systems and value not in candidate_systems:
                            candidate_systems.append(value)
            if not candidate_systems and len(systems) == 1:
                candidate_systems = list(systems)
            if candidate_systems:
                selector = candidate_systems[0] if len(candidate_systems) == 1 else candidate_systems
                result = sdata.filter_by_coordinate_system(selector)
                if result is None:
                    raise ValueError(f"coordinate-system split for {values!r} returned no data")
                result = _copy_sdata_with_table(result, table_key, subset_table)
            else:
                logger.warning(
                    "No coordinate system maps unambiguously to %s=%s; retaining spatial elements and filtering instances",
                    split_by,
                    ",".join(values),
                )
                result = _sdata_for_table_subset(sdata, table_key, subset_table, args.subset_mode)
            if individual:
                filename = f"{safe_name(values[0])}.zarr"
            else:
                filename = f"{safe_name(split_by)}_selected_{safe_name('_'.join(values))}.zarr"
            output_path = os.path.join(args.output_dir, filename)
            result.write(output_path, overwrite=True)
            export_annotation_csvs(
                subset_table, output_path, split_by, args.annotation_format,
                args.annotation_cols, sdata=result,
            )
        return

    if "group" not in adata.obs.columns:
        raise ValueError("group column not found in table obs")
    region_key, _ = _table_keys(adata)
    if not region_key or region_key not in adata.obs.columns:
        logger.warning("No valid region_key is available; group split will subset the table only")
        one_to_one = False
    else:
        mapping_counts = adata.obs.groupby(region_key, observed=True)["group"].nunique(dropna=False)
        one_to_one = bool((mapping_counts <= 1).all())
    obs_values = adata.obs["group"].astype(str)
    groups = parse_barcode_groups(args.barcodes)
    individual = not groups
    if individual:
        groups = [[str(value)] for value in pd.unique(obs_values)]
    available = set(pd.unique(obs_values))
    for group_values in groups:
        missing = [value for value in group_values if value not in available]
        if missing:
            raise ValueError(f"values not found in group: {', '.join(missing)}")
        mask = obs_values.isin(group_values)
        subset_table = prepare_subset_table(adata, mask.to_numpy())
        if individual:
            filename = f"group_{safe_name(group_values[0])}.zarr"
        else:
            filename = f"group_selected_{safe_name('_'.join(group_values))}.zarr"
        output_path = os.path.join(args.output_dir, filename)
        regions = (
            pd.unique(subset_table.obs[region_key].astype(str)).tolist()
            if region_key and region_key in subset_table.obs.columns else []
        )
        if one_to_one and all(region in sdata.coordinate_systems for region in regions):
            pieces = [sdata.filter_by_coordinate_system(region) for region in regions]
            pieces = [piece for piece in pieces if piece is not None]
            if not pieces:
                raise ValueError(f"no coordinate systems found for group {group_values!r}")
            result = pieces[0] if len(pieces) == 1 else spd.concatenate(pieces, concatenate_tables=True)
            # Make the selected table exact even if a coordinate system contains
            # observations not represented by the group column.
            result = _copy_sdata_with_table(result, table_key, subset_table)
            result.write(output_path, overwrite=True)
            export_annotation_csvs(subset_table, output_path, "group", args.annotation_format,
                                   args.annotation_cols, sdata=result)
        else:
            if region_key and region_key in adata.obs.columns and not one_to_one:
                logger.warning("At least one region contains multiple groups; filtering table instances exactly")
            _write_subset(subset_table, output_path, input_kind, sdata, table_key, args, "group")


def _csv_header_row(path):
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as handle:
        for number, line in enumerate(handle):
            fields = {_normalise_column(field) for field in line.rstrip("\r\n").split(",")}
            if "barcode" in fields or "cellid" in fields or {"x", "y"}.issubset(fields):
                return number
    return 0


def _read_roi_csv(path):
    return pd.read_csv(path, skiprows=_csv_header_row(path), sep=",", engine="python", on_bad_lines="skip")


def _selection_name_from_preamble(path):
    header_row = _csv_header_row(path)
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as handle:
        for number, line in enumerate(handle):
            if number >= header_row:
                break
            match = re.match(r"\s*#\s*Selection\s+name\s*:\s*(.+?)\s*$", line, flags=re.IGNORECASE)
            if match:
                return match.group(1).strip()
    return ""


def _find_column(df, requested):
    normal = {_normalise_column(column): column for column in df.columns}
    return normal.get(_normalise_column(requested))


def load_roi_definitions(roi_csv, roi_label_col=""):
    if not roi_csv:
        raise ValueError("--roi_csv is required for ROI splitting")
    if os.path.isdir(roi_csv):
        paths = sorted(
            os.path.join(roi_csv, filename)
            for filename in os.listdir(roi_csv)
            if filename.lower().endswith(".csv")
        )
    else:
        paths = [roi_csv]
    if not paths:
        raise ValueError("no CSV files found for ROI splitting")

    definitions = []
    for path in paths:
        df = _read_roi_csv(path)
        if df.empty:
            raise ValueError(f"ROI CSV {path!r} contains no data rows")
        columns = {_normalise_column(column): column for column in df.columns}
        x_col, y_col = columns.get("x"), columns.get("y")
        id_candidate = columns.get("barcode") or columns.get("cellid")
        stem = Path(path).stem
        preamble_name = _selection_name_from_preamble(path)
        if id_candidate is None and x_col is not None and y_col is not None:
            label_col = _find_column(df, roi_label_col) if roi_label_col else columns.get("selection")
            if roi_label_col and label_col is None:
                raise ValueError(f"roi_label_col {roi_label_col!r} not found in {path}")
            groups = [(preamble_name or stem, df)] if label_col is None else df.groupby(label_col, sort=False, dropna=False)
            for roi_name, roi_df in groups:
                coords = roi_df[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna().to_numpy()
                if len({tuple(point) for point in coords}) < 3:
                    raise ValueError(f"polygon ROI {roi_name!r} must contain at least three distinct vertices")
                if tuple(coords[0]) != tuple(coords[-1]):
                    raise ValueError(f"polygon ROI {roi_name!r} is not closed (first and last coordinates differ)")
                polygon = Polygon(coords)
                if polygon.is_empty or polygon.area <= 0 or not polygon.is_valid:
                    raise ValueError(f"polygon ROI {roi_name!r} is not a valid closed polygon")
                definitions.append({"kind": "polygon", "name": str(roi_name), "polygon": polygon,
                                    "hint": "xenium"})
            continue

        id_col = id_candidate
        if id_col is None:
            # A one-column barcode list is a common Loupe export.
            if len(df.columns) == 1:
                id_col = df.columns[0]
            else:
                raise ValueError(f"barcode or Cell ID column not found in {path}")
        hint = "loupe" if _normalise_column(id_col) == "barcode" else "xenium"
        excluded = {id_col}
        scope_col = None
        for candidate in ("source_region", "sample_id", "sample"):
            found = _find_column(df, candidate)
            if found is not None:
                scope_col = found
                excluded.add(found)
                break
        if roi_label_col:
            label_col = _find_column(df, roi_label_col)
            if label_col is None:
                raise ValueError(f"roi_label_col {roi_label_col!r} not found in {path}")
        else:
            preferred = ("roi", "group", "selection", "category", "grouped_annotation", "annotation_name", "region")
            label_col = next((_find_column(df, value) for value in preferred if _find_column(df, value) not in excluded), None)
            statistic_prefixes = (
                "cluster", "transcript", "totaltranscript", "area", "cellarea",
                "count", "density",
            )
            remaining = [
                column for column in df.columns
                if column not in excluded
                and _normalise_column(column) not in {"x", "y"}
                and not _normalise_column(column).startswith(statistic_prefixes)
            ]
            if label_col is None and len(remaining) == 1:
                label_col = remaining[0]
            elif label_col is None and len(remaining) > 1:
                raise ValueError(f"multiple ROI label columns found in {path}; specify --roi_label_col")
        groups = [(preamble_name or stem, df)] if label_col is None else df.groupby(label_col, sort=False, dropna=False)
        for roi_name, roi_df in groups:
            values = roi_df[id_col].dropna().astype(str).str.strip()
            values = values[values.ne("")]
            if values.empty:
                raise ValueError(f"ROI {roi_name!r} contains no cell/barcode IDs")
            definition = {"kind": "ids", "name": str(roi_name), "ids": values.tolist(), "hint": hint}
            if scope_col is not None:
                definition["regions"] = roi_df.loc[values.index, scope_col].astype(str).tolist()
            definitions.append(definition)
    return definitions


def _obs_identifier(adata):
    _, instance_key = _table_keys(adata)
    normal = {_normalise_column(column): column for column in adata.obs.columns}
    if instance_key and instance_key in adata.obs.columns:
        return adata.obs[instance_key].astype(str), instance_key
    if "cellid" in normal:
        column = normal["cellid"]
        return adata.obs[column].astype(str), column
    return pd.Series(adata.obs_names.astype(str), index=adata.obs_names), adata.obs_names.name or "obs_names"


def _mask_for_id_roi(adata, definition, roi_region=""):
    ids, identifier_name = _obs_identifier(adata)
    requested = pd.Index(definition["ids"]).astype(str)
    region_key, _ = _table_keys(adata)
    region_series = adata.obs[region_key].astype(str) if region_key and region_key in adata.obs.columns else None
    if roi_region:
        if region_series is None:
            raise ValueError("--roi_region was provided but the table has no valid region_key")
        if roi_region not in set(region_series):
            raise ValueError(f"roi_region {roi_region!r} not found in {region_key}")
        mask = ids.isin(requested) & region_series.eq(roi_region)
    elif "regions" in definition:
        if region_series is None:
            raise ValueError("ROI CSV contains sample/region information but the table has no valid region_key")
        pairs = set(zip(requested, map(str, definition["regions"])))
        mask = pd.Series([(identifier, region) in pairs for identifier, region in zip(ids, region_series)], index=adata.obs_names)
    else:
        mask = ids.isin(requested)
        duplicated = ids[mask].duplicated(keep=False)
        if duplicated.any() and region_series is not None and region_series[mask].nunique() > 1:
            raise ValueError(
                f"ROI IDs are duplicated across integrated regions for {identifier_name}; "
                "provide --roi_region or sample information in the CSV"
            )
    matched_ids = set(ids[mask])
    unmatched = [value for value in pd.unique(requested) if value not in matched_ids]
    logger.info("ROI %s: requested=%d, matched=%d, unmatched=%d", definition["name"],
                len(pd.unique(requested)), int(mask.sum()), len(unmatched))
    if not mask.any():
        raise ValueError(f"ROI {definition['name']!r} matched no observations")
    if unmatched:
        logger.warning("ROI %s contains %d unmatched IDs (first values: %s)", definition["name"],
                       len(unmatched), ", ".join(unmatched[:5]))
    return mask.to_numpy()


def _polygon_mask_from_obsm(adata, polygon):
    if "spatial" not in adata.obsm:
        return None
    coords = adata.obsm["spatial"]
    if coords.shape[1] < 2:
        return None
    return pd.Series([polygon.covers(Point(float(x), float(y))) for x, y in coords[:, :2]], index=adata.obs_names)


def _polygon_mask(adata, sdata, polygon, coordinate_system):
    if sdata is None:
        mask = _polygon_mask_from_obsm(adata, polygon)
        if mask is None:
            raise ValueError("polygon ROI requires SpatialData elements or adata.obsm['spatial']")
        return mask.to_numpy()

    region_key, instance_key = _table_keys(adata)
    if not instance_key or instance_key not in adata.obs.columns:
        raise ValueError("polygon ROI requires a valid spatialdata_attrs.instance_key")
    selected_pairs = set()
    annotated_regions = _spatialdata_attrs(adata).get("region", [])
    if isinstance(annotated_regions, str):
        annotated_regions = [annotated_regions]
    for element_name in annotated_regions:
        if element_name in sdata.shapes:
            element = spd.transform(sdata.shapes[element_name], to_coordinate_system=coordinate_system)
            geometries = element.geometry
            if "radius" in element.columns:
                geometries = geometries.buffer(pd.to_numeric(element["radius"], errors="coerce").fillna(0))
            selected = pd.Series(
                [polygon.covers(geometry) for geometry in geometries],
                index=element.index,
            )
            selected_pairs.update((str(element_name), str(index)) for index in element.index[selected])
        elif element_name in sdata.points:
            element = spd.transform(sdata.points[element_name], to_coordinate_system=coordinate_system)
            frame = element.compute() if hasattr(element, "compute") else element
            selected = [polygon.covers(Point(float(x), float(y))) for x, y in zip(frame["x"], frame["y"])]
            selected_pairs.update((str(element_name), str(index)) for index in frame.index[selected])
    if not selected_pairs:
        raise ValueError("polygon ROI matched no Shapes or Points instances")
    instances = adata.obs[instance_key].astype(str)
    if region_key and region_key in adata.obs.columns:
        regions = adata.obs[region_key].astype(str)
        mask = pd.Series([(region, instance) in selected_pairs for region, instance in zip(regions, instances)],
                         index=adata.obs_names)
    else:
        selected_ids = {instance for _, instance in selected_pairs}
        mask = instances.isin(selected_ids)
    if not mask.any():
        raise ValueError("polygon ROI spatial instances did not match table observations")
    return mask.to_numpy()


def split_roi(adata, input_kind, sdata, table_key, args):
    definitions = load_roi_definitions(args.roi_csv, args.roi_label_col)
    coordinate_system = None
    if any(definition["kind"] == "polygon" for definition in definitions):
        if input_kind == "zarr":
            coordinate_system = resolve_coordinate_system(
                sdata, args.coordinate_system, args.shape_elements, args.roi_region
            )
    for definition in definitions:
        if definition["kind"] == "ids":
            mask = _mask_for_id_roi(adata, definition, args.roi_region)
            polygon = None
        else:
            mask = _polygon_mask(adata, sdata, definition["polygon"], coordinate_system)
            polygon = definition["polygon"]
            logger.info("ROI %s: matched=%d observations by polygon", definition["name"], int(mask.sum()))
        subset = prepare_subset_table(adata, mask)
        output_path = os.path.join(args.output_dir, f"ROI_{safe_name(definition['name'])}.{input_kind}")
        _write_subset(
            subset, output_path, input_kind, sdata, table_key, args, "ROI",
            roi_name=definition["name"], roi_hint=definition["hint"], polygon=polygon,
            coordinate_system=coordinate_system,
        )
        logger.info("Wrote ROI %s (%d observations) to %s", definition["name"], subset.n_obs, output_path)


def crop_image(sdata, table_key, args):
    if not (args.min_x < args.max_x and args.min_y < args.max_y):
        raise ValueError("image bounds must satisfy min_x < max_x and min_y < max_y")
    coordinate_system = resolve_coordinate_system(sdata, args.coordinate_system, args.shape_elements)
    subset = spd.bounding_box_query(
        sdata,
        min_coordinate=[args.min_x, args.min_y],
        max_coordinate=[args.max_x, args.max_y],
        axes=("x", "y"),
        target_coordinate_system=coordinate_system,
    )
    if subset is None:
        raise ValueError("image crop did not intersect the SpatialData object")
    image_id = f"{args.min_x:g}_{args.max_x:g}_{args.min_y:g}_{args.max_y:g}"
    output_path = os.path.join(args.output_dir, f"spatial{image_id}.zarr")
    subset.write(output_path, overwrite=True)
    if table_key in subset.tables:
        export_annotation_csvs(
            subset.tables[table_key], output_path, "image", args.annotation_format,
            args.annotation_cols, sdata=subset, roi_name=f"spatial_{image_id}",
        )
    _save_crop_preview(subset, coordinate_system, os.path.join(args.output_dir, f"{image_id}_shape.png"))


def _save_crop_preview(sdata, coordinate_system, output_path):
    """Best-effort preview supporting images, shapes, and points-only objects."""
    try:
        render = sdata.pl
        rendered = False
        if len(sdata.images):
            render = render.render_images()
            rendered = True
        if len(sdata.shapes):
            render = render.render_shapes()
            rendered = True
        if len(sdata.points):
            render = render.render_points()
            rendered = True
        if not rendered:
            logger.warning("Crop is valid but contains no plottable images, shapes, or points")
            return
        _, ax = plt.subplots(1, 1, figsize=(12, 10))
        render.pl.show(ax=ax, coordinate_systems=coordinate_system)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as exc:
        plt.close("all")
        logger.warning("Crop was written, but preview rendering failed: %s", exc)


def main(argv=None):
    args = build_parser().parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    log_step(logger, 1, 4, f"preparing split by {args.split_by}")
    input_kind = _detect_input(args.INPUT_FILE)
    log_step(logger, 2, 4, f"loading input file {args.INPUT_FILE}")
    if input_kind == "h5ad":
        sdata = None
        table_key = None
        adata = sc.read_h5ad(args.INPUT_FILE)
    else:
        sdata = spd.read_zarr(args.INPUT_FILE)
        table_key = _select_table(sdata, args.table_key)
        adata = sdata.tables[table_key]
    logger.info("Loaded %d observations for splitting", adata.n_obs)

    split_lower = args.split_by.lower()
    log_step(logger, 3, 4, f"running {args.split_by} split")
    if split_lower in {"image", "images"}:
        if input_kind != "zarr":
            raise ValueError("image splitting requires a SpatialData Zarr input")
        crop_image(sdata, table_key, args)
    elif split_lower in {"roi", "rois"}:
        split_roi(adata, input_kind, sdata, table_key, args)
    elif args.split_by in adata.obs.columns:
        if split_lower in {"sample", "samples", "region", "group"}:
            split_sample_or_group(adata, input_kind, sdata, table_key, args)
        else:
            split_metadata(adata, input_kind, sdata, table_key, args)
    else:
        raise ValueError(f"unknown split_by {args.split_by!r}; use ROI, image, or an obs column")
    log_step(logger, 4, 4, f"split output saved to {args.output_dir}")


if __name__ == "__main__":
    main()
