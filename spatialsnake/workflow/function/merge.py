"""Merge independent SpatialData samples or overlay subset annotations.

This module is intentionally executed directly by ``spatialsnake useful_tool``.
It does not perform batch correction, spatial registration, or expression
normalization.
"""

from __future__ import annotations

import argparse
import copy
import os
import shutil
import uuid
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import spatialdata as sd
from pandas.api.types import (
    is_bool_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_string_dtype,
)

from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("useful_merge")

AUTO_VALUES = {"", "auto", "none", "null"}
CSV_CELL_CANDIDATES = ("barcode", "cell_id", "cellid", "cell_barcode")
CSV_REGION_CANDIDATES = ("region", "sample", "sample_id", "library_id")
CSV_ANNOTATION_CANDIDATES = (
    "grouped_annotation",
    "celltype",
    "annotation",
    "cluster_id",
    "group",
)
ZARR_ANNOTATION_CANDIDATES = ("sub_celltype", "celltype", "clusters")


def _is_auto(value: Any) -> bool:
    return str(value).strip().lower() in AUTO_VALUES


def _split_values(value: Any) -> list[str]:
    if _is_auto(value):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _read_sdata(path: str) -> sd.SpatialData:
    if not os.path.isdir(path):
        raise FileNotFoundError(f"SpatialData Zarr directory does not exist: {path}")
    return sd.read_zarr(path)


def _copy_sdata_container(source: sd.SpatialData) -> sd.SpatialData:
    """Copy the container mappings while keeping lazy spatial elements lazy."""

    return sd.SpatialData(
        images=dict(source.images),
        labels=dict(source.labels),
        points=dict(source.points),
        shapes=dict(source.shapes),
        tables=dict(source.tables),
        attrs=copy.deepcopy(source.attrs),
    )


def _resolve_table_key(
    sdata: sd.SpatialData,
    requested: str,
    source: str,
    *,
    allow_single_fallback: bool = False,
) -> str:
    keys = list(sdata.tables)
    if not keys:
        raise ValueError(f"No tables were found in {source}.")
    if not _is_auto(requested):
        if requested in sdata.tables:
            return requested
        if allow_single_fallback and len(keys) == 1:
            return keys[0]
        raise ValueError(
            f"Table {requested!r} was not found in {source}; available tables: {keys}."
        )
    if len(keys) != 1:
        raise ValueError(
            f"{source} contains multiple tables {keys}; set --table_key explicitly."
        )
    return keys[0]


def _spatial_attrs(table: Any, source: str) -> dict[str, Any]:
    attrs = table.uns.get("spatialdata_attrs")
    if not isinstance(attrs, dict):
        raise ValueError(f"Table in {source} has no valid spatialdata_attrs metadata.")
    region_key = attrs.get("region_key")
    instance_key = attrs.get("instance_key")
    if not region_key or region_key not in table.obs:
        raise ValueError(
            f"Table in {source} has an invalid region_key {region_key!r}."
        )
    if not instance_key or instance_key not in table.obs:
        raise ValueError(
            f"Table in {source} has an invalid instance_key {instance_key!r}."
        )
    return attrs


def _validate_table_links(sdata: sd.SpatialData, table_key: str, source: str) -> None:
    table = sdata.tables[table_key]
    attrs = _spatial_attrs(table, source)
    region_key = attrs["region_key"]
    instance_key = attrs["instance_key"]
    if table.obs[region_key].isna().any() or table.obs[instance_key].isna().any():
        raise ValueError(f"Table {table_key!r} in {source} has missing spatial link values.")
    spatial_names = set(sdata.shapes) | set(sdata.points) | set(sdata.labels)
    table_regions = set(table.obs[region_key].astype(str))
    missing = sorted(table_regions.difference(spatial_names))
    if missing:
        raise ValueError(
            f"Table {table_key!r} in {source} references absent spatial element(s): {missing}."
        )


def _dtype_family(dtype: Any) -> str:
    if isinstance(dtype, pd.CategoricalDtype) or is_string_dtype(dtype) or dtype == object:
        return "text"
    if is_bool_dtype(dtype):
        return "bool"
    if is_numeric_dtype(dtype):
        return "numeric"
    if is_datetime64_any_dtype(dtype):
        return "datetime"
    return str(dtype)


def _validate_obs_dtypes(tables: Sequence[Any]) -> None:
    common = set(tables[0].obs.columns)
    for table in tables[1:]:
        common.intersection_update(table.obs.columns)
    for column in sorted(common):
        families = {_dtype_family(table.obs[column].dtype) for table in tables}
        if len(families) > 1:
            raise ValueError(
                f"Observation column {column!r} has incompatible dtype families: {sorted(families)}."
            )


def _validate_layers_and_raw(tables: Sequence[Any]) -> None:
    layer_sets = [set(table.layers) for table in tables]
    if any(keys != layer_sets[0] for keys in layer_sets[1:]):
        raise ValueError(
            "Input tables contain different layer keys; harmonize layers before sample merge. "
            f"Got: {[sorted(keys) for keys in layer_sets]}."
        )
    raw_presence = [table.raw is not None for table in tables]
    if len(set(raw_presence)) > 1:
        raise ValueError("Either all sample tables must contain .raw or none of them may contain .raw.")


def _align_obsm(tables: Sequence[Any]) -> None:
    """Align named embeddings before AnnData's single join strategy is applied."""

    all_keys: set[str] = set()
    for table in tables:
        all_keys.update(table.obsm.keys())

    for key in sorted(all_keys):
        values = [table.obsm[key] if key in table.obsm else None for table in tables]
        present = [value for value in values if value is not None]
        frames = [isinstance(value, pd.DataFrame) for value in present]
        if any(frames):
            if not all(frames):
                raise ValueError(
                    f"obsm[{key!r}] mixes DataFrame and unnamed array representations."
                )
            union: list[str] = []
            for frame in present:
                for column in frame.columns.astype(str):
                    if column not in union:
                        union.append(column)
            for table, value in zip(tables, values):
                if value is None:
                    table.obsm[key] = pd.DataFrame(
                        np.nan,
                        index=table.obs_names.copy(),
                        columns=union,
                    )
                else:
                    frame = value.copy()
                    frame.columns = frame.columns.astype(str)
                    table.obsm[key] = frame.reindex(columns=union)
            continue

        if len(present) != len(tables):
            raise ValueError(
                f"Unnamed obsm[{key!r}] is missing from one or more samples and cannot be safely aligned."
            )
        trailing_shapes = {tuple(value.shape[1:]) for value in present}
        if len(trailing_shapes) != 1:
            raise ValueError(
                f"Unnamed obsm[{key!r}] has incompatible dimensions: {sorted(trailing_shapes)}."
            )


def _prepare_obs_names(tables: Sequence[Any], sample_ids: Sequence[str]) -> bool:
    """Add sample suffixes only when observation IDs collide across inputs."""

    all_names = [
        str(name)
        for table in tables
        for name in table.obs_names
    ]
    if len(all_names) == len(set(all_names)):
        logger.info("Observation IDs are already unique across inputs; preserving obs_names.")
        return False

    for table, sample_id in zip(tables, sample_ids):
        suffix = f"-{sample_id}"
        table.obs_names = pd.Index(
            [
                name if str(name).endswith(suffix) else f"{name}{suffix}"
                for name in table.obs_names.astype(str)
            ],
            name=table.obs_names.name,
        )
        # Keep DataFrame-valued obsm indices synchronized with AnnData.obs_names.
        for key in table.obsm:
            if isinstance(table.obsm[key], pd.DataFrame):
                frame = table.obsm[key].copy()
                frame.index = table.obs_names.copy()
                table.obsm[key] = frame

    renamed = [
        str(name)
        for table in tables
        for name in table.obs_names
    ]
    if len(renamed) != len(set(renamed)):
        raise ValueError(
            "Observation IDs remain duplicated after adding sample suffixes. "
            "Use distinct --sample_ids or provide unique input observation IDs."
        )
    logger.info("Duplicated observation IDs detected; added sample suffixes once.")
    return True


def _spatial_element_names_are_unique(sdatas: Sequence[sd.SpatialData]) -> bool:
    names = [
        name
        for sdata in sdatas
        for _, name, _ in sdata.gen_spatial_elements()
    ]
    return len(names) == len(set(names))


def _resolve_sample_ids(
    paths: Sequence[str],
    tables: Sequence[Any],
    explicit: str,
    sample_col: str,
) -> list[str]:
    requested = _split_values(explicit)
    if requested and len(requested) != len(paths):
        raise ValueError(
            f"--sample_ids supplied {len(requested)} ID(s) for {len(paths)} input object(s)."
        )
    if requested:
        sample_ids = requested
    else:
        sample_ids = []
        for path, table in zip(paths, tables):
            sample_id = ""
            if sample_col in table.obs:
                values = table.obs[sample_col].dropna().astype(str)
                values = values[values.str.strip().ne("")].unique().tolist()
                if len(values) == 1:
                    sample_id = values[0]
            if not sample_id:
                name = Path(path).name
                sample_id = name[:-5] if name.lower().endswith(".zarr") else name
            sample_ids.append(sample_id)
    if any(not str(value).strip() for value in sample_ids):
        raise ValueError("Sample IDs must be non-empty.")
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError(f"Sample IDs must be unique; got {sample_ids}.")
    return sample_ids


def _apply_label_policy(
    tables: Sequence[Any],
    sample_ids: Sequence[str],
    label_cols: str,
    policy: str,
) -> None:
    columns = _split_values(label_cols)
    if policy == "preserve":
        return
    if not columns:
        raise ValueError(f"label_policy={policy!r} requires --label_cols.")

    counters = {column: 0 for column in columns}
    for table, sample_id in zip(tables, sample_ids):
        for column in columns:
            if column not in table.obs:
                raise ValueError(f"Label column {column!r} is missing from sample {sample_id!r}.")
            values = table.obs[column]
            if policy == "prefix":
                mapped = values.astype(object).copy()
                mask = values.notna()
                mapped.loc[mask] = sample_id + ":" + values.loc[mask].astype(str)
                table.obs[column] = pd.Categorical(mapped)
            elif policy == "offset":
                categories = pd.unique(values.dropna())
                mapping = {
                    value: index + counters[column]
                    for index, value in enumerate(categories)
                }
                table.obs[column] = values.map(mapping).astype("Int64")
                counters[column] += len(categories)


def _feature_counts(tables: Sequence[Any]) -> tuple[int, int]:
    intersection = set(tables[0].var_names.astype(str))
    union = set(intersection)
    for table in tables[1:]:
        names = set(table.var_names.astype(str))
        intersection.intersection_update(names)
        union.update(names)
    return len(intersection), len(union)


def _append_merge_history(uns: Any, entry: dict[str, Any]) -> None:
    """Store history as a Zarr-serializable mapping, not a list of mappings."""

    existing = uns.get("merge_history", {})
    if isinstance(existing, dict):
        history = copy.deepcopy(existing)
    elif isinstance(existing, list):
        history = {
            f"run_{index + 1:04d}": copy.deepcopy(value)
            for index, value in enumerate(existing)
        }
    else:
        history = {"run_0001_previous": str(existing)}
    history[f"run_{len(history) + 1:04d}"] = entry
    uns["merge_history"] = history


def _sample_report_row(
    path: str,
    sample_id: str,
    sdata: sd.SpatialData,
    table_key: str,
    intersection: int,
    union: int,
) -> dict[str, Any]:
    table = sdata.tables[table_key]
    attrs = _spatial_attrs(table, path)
    regions = sorted(table.obs[attrs["region_key"]].astype(str).unique())
    return {
        "mode": "sample",
        "source": path,
        "sample_id": sample_id,
        "table_key": table_key,
        "n_obs": table.n_obs,
        "n_vars": table.n_vars,
        "regions": ";".join(regions),
        "images": len(sdata.images),
        "labels": len(sdata.labels),
        "points": len(sdata.points),
        "shapes": len(sdata.shapes),
        "obs_fields": ";".join(map(str, table.obs.columns)),
        "obsm_fields": ";".join(map(str, table.obsm.keys())),
        "feature_intersection": intersection,
        "feature_union": union,
    }


def merge_samples(
    input_paths: Sequence[str],
    *,
    table_key: str = "",
    sample_ids: str = "",
    sample_col: str = "sample",
    feature_join: str = "inner",
    label_cols: str = "",
    label_policy: str = "preserve",
) -> tuple[sd.SpatialData, str, list[dict[str, Any]]]:
    if len(input_paths) < 2:
        raise ValueError("merge_by=sample requires at least two SpatialData inputs.")

    source_sdatas = [_read_sdata(path) for path in input_paths]
    selected_keys = [
        _resolve_table_key(sdata, table_key, path)
        for path, sdata in zip(input_paths, source_sdatas)
    ]
    for path, sdata, key in zip(input_paths, source_sdatas, selected_keys):
        _validate_table_links(sdata, key, path)

    tables = [
        sdata.tables[key].copy()
        for sdata, key in zip(source_sdatas, selected_keys)
    ]
    resolved_ids = _resolve_sample_ids(input_paths, tables, sample_ids, sample_col)
    _validate_obs_dtypes(tables)
    _validate_layers_and_raw(tables)
    _align_obsm(tables)
    _apply_label_policy(tables, resolved_ids, label_cols, label_policy)
    obs_names_renamed = _prepare_obs_names(tables, resolved_ids)

    intersection, union = _feature_counts(tables)
    if feature_join == "inner" and intersection == 0:
        raise ValueError("Input samples have no shared genes; inner feature join would be empty.")
    logger.info(
        "Feature alignment: intersection=%d, union=%d, strategy=%s",
        intersection,
        union,
        feature_join,
    )

    merged_table_key = table_key if not _is_auto(table_key) else selected_keys[0]
    merge_sources: dict[str, Any] = {}
    prepared: dict[str, sd.SpatialData] = {}
    report: list[dict[str, Any]] = []

    for path, sample_id, source, selected_key, table in zip(
        input_paths, resolved_ids, source_sdatas, selected_keys, tables
    ):
        if sample_col not in table.obs:
            table.obs[sample_col] = sample_id
        else:
            values = table.obs[sample_col].astype(object).copy()
            missing = table.obs[sample_col].isna() | table.obs[sample_col].astype(str).str.strip().eq("")
            values.loc[missing] = sample_id
            table.obs[sample_col] = values
        table.obs["_merge_source"] = sample_id

        if not table.obs_names.is_unique:
            raise ValueError(f"Observation IDs are not unique within sample {sample_id!r}.")

        non_spatial_uns = {
            key: copy.deepcopy(value)
            for key, value in table.uns.items()
            if key != "spatialdata_attrs"
        }
        merge_sources[sample_id] = non_spatial_uns

        working = _copy_sdata_container(source)
        working.tables[selected_key] = table
        if selected_key != merged_table_key:
            if merged_table_key in working.tables:
                raise ValueError(
                    f"Cannot rename selected table {selected_key!r} to {merged_table_key!r} in {path}; "
                    "that table name already exists."
                )
            del working.tables[selected_key]
            working.tables[merged_table_key] = table

        for other_key in list(working.tables):
            if other_key == merged_table_key:
                continue
            renamed = f"{other_key}__{sample_id}"
            if renamed in working.tables:
                raise ValueError(f"Table rename collision in {path}: {renamed!r}.")
            other_table = working.tables[other_key]
            del working.tables[other_key]
            working.tables[renamed] = other_table

        prepared[sample_id] = working
        report.append(
            _sample_report_row(path, sample_id, source, selected_key, intersection, union)
        )

    spatial_names_unique = _spatial_element_names_are_unique(list(prepared.values()))
    concatenate_inputs: Sequence[sd.SpatialData] | dict[str, sd.SpatialData]
    if spatial_names_unique:
        logger.info("Spatial element names are already unique; preserving element and region names.")
        concatenate_inputs = list(prepared.values())
    else:
        logger.info("Duplicated spatial element names detected; adding sample suffixes once.")
        concatenate_inputs = prepared

    merged = sd.concatenate(
        concatenate_inputs,
        concatenate_tables=True,
        obs_names_make_unique=False,
        join=feature_join,
        merge="same",
        uns_merge="same",
        pairwise=True,
    )
    if merged_table_key not in merged.tables:
        raise RuntimeError(f"Merged table {merged_table_key!r} was not created.")
    merged_table = merged.tables[merged_table_key]
    expected_obs = sum(table.n_obs for table in tables)
    if merged_table.n_obs != expected_obs:
        raise RuntimeError(
            f"Merged table contains {merged_table.n_obs} observations; expected {expected_obs}."
        )
    merged_table.uns["merge_sources"] = merge_sources
    _append_merge_history(
        merged_table.uns,
        {
            "mode": "sample",
            "sources": list(resolved_ids),
            "feature_join": feature_join,
            "feature_intersection": intersection,
            "feature_union": union,
            "label_policy": label_policy,
            "label_cols": _split_values(label_cols),
            "obs_names_renamed": obs_names_renamed,
            "spatial_names_renamed": not spatial_names_unique,
        },
    )
    merged.tables[merged_table_key] = merged_table
    _validate_table_links(merged, merged_table_key, "merged output")
    return merged, merged_table_key, report


def _parse_csv_inputs(annotation_csv: str) -> list[str]:
    paths: list[str] = []
    for part in _split_values(annotation_csv):
        if os.path.isdir(part):
            paths.extend(
                str(path)
                for path in sorted(Path(part).iterdir())
                if path.suffix.lower() == ".csv"
            )
        else:
            paths.append(part)
    return paths


def _choose_column(
    columns: Iterable[Any],
    requested: str,
    candidates: Sequence[str],
    source: str,
    purpose: str,
) -> Any:
    columns = list(columns)
    lower = {str(column).strip().lower(): column for column in columns}
    if not _is_auto(requested):
        if requested in columns:
            return requested
        match = lower.get(str(requested).strip().lower())
        if match is not None:
            return match
        raise ValueError(
            f"{purpose} column {requested!r} was not found in {source}; columns={list(map(str, columns))}."
        )
    matches = [lower[candidate] for candidate in candidates if candidate in lower]
    matches = list(dict.fromkeys(matches))
    if len(matches) != 1:
        raise ValueError(
            f"Could not unambiguously infer the {purpose} column in {source}; "
            f"candidates found={list(map(str, matches))}. Set the corresponding parameter explicitly."
        )
    return matches[0]


def _normalized_text(series: pd.Series, source: str, purpose: str) -> pd.Series:
    if series.isna().any():
        raise ValueError(f"{source} contains missing {purpose} values.")
    normalized = series.astype(str).str.strip()
    if normalized.eq("").any():
        raise ValueError(f"{source} contains empty {purpose} values.")
    return normalized


def _base_match_context(
    table: Any,
    source: str,
    input_cell_col: str,
    input_region_col: str,
) -> dict[str, Any]:
    attrs = _spatial_attrs(table, source)
    cell_col = attrs["instance_key"] if _is_auto(input_cell_col) else input_cell_col
    region_col = attrs["region_key"] if _is_auto(input_region_col) else input_region_col
    if cell_col not in table.obs:
        raise ValueError(f"Base cell ID column {cell_col!r} was not found in {source}.")
    if region_col not in table.obs:
        raise ValueError(f"Base region column {region_col!r} was not found in {source}.")

    cells = _normalized_text(table.obs[cell_col], source, "cell ID")
    regions = _normalized_text(table.obs[region_col], source, "region")
    use_region = regions.nunique(dropna=False) > 1 or cells.duplicated().any()
    keys = list(zip(regions, cells)) if use_region else cells.tolist()
    if len(set(keys)) != len(keys):
        raise ValueError(
            f"Base matching keys are not unique in {source}; use valid region and instance columns."
        )
    return {
        "cell_col": cell_col,
        "region_col": region_col,
        "cells": cells,
        "regions": regions,
        "keys": keys,
        "use_region": use_region,
    }


def _records_from_csv(
    path: str,
    *,
    annotation_col: str,
    csv_cell_col: str,
    csv_region_col: str,
    use_region: bool,
) -> tuple[list[Any], list[Any], int]:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Annotation CSV does not exist: {path}")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"Annotation CSV is empty: {path}")
    cell_col = _choose_column(
        frame.columns, csv_cell_col, CSV_CELL_CANDIDATES, path, "cell ID"
    )
    label_col = _choose_column(
        frame.columns,
        annotation_col,
        CSV_ANNOTATION_CANDIDATES,
        path,
        "annotation",
    )
    cells = _normalized_text(frame[cell_col], path, "cell ID")
    if use_region:
        region_col = _choose_column(
            frame.columns,
            csv_region_col,
            CSV_REGION_CANDIDATES,
            path,
            "region",
        )
        regions = _normalized_text(frame[region_col], path, "region")
        keys = list(zip(regions, cells))
    else:
        keys = cells.tolist()
    valid = frame[label_col].notna() & frame[label_col].astype(str).str.strip().ne("")
    invalid_labels = int((~valid).sum())
    valid_array = valid.to_numpy()
    selected_keys = [key_value for key_value, keep in zip(keys, valid_array) if keep]
    return selected_keys, frame.loc[valid, label_col].tolist(), invalid_labels


def _records_from_zarr(
    path: str,
    *,
    table_key: str,
    annotation_col: str,
    use_region: bool,
) -> tuple[list[Any], list[Any], int]:
    source = _read_sdata(path)
    key = _resolve_table_key(
        source,
        table_key,
        path,
        allow_single_fallback=True,
    )
    table = source.tables[key]
    attrs = _spatial_attrs(table, path)
    label_col = _choose_column(
        table.obs.columns,
        annotation_col,
        ZARR_ANNOTATION_CANDIDATES,
        path,
        "annotation",
    )
    cells = _normalized_text(table.obs[attrs["instance_key"]], path, "cell ID")
    if use_region:
        regions = _normalized_text(table.obs[attrs["region_key"]], path, "region")
        keys = list(zip(regions, cells))
    else:
        keys = cells.tolist()
    labels = table.obs[label_col]
    valid = labels.notna() & labels.astype(str).str.strip().ne("")
    invalid_labels = int((~valid).sum())
    valid_array = valid.to_numpy()
    selected_keys = [key_value for key_value, keep in zip(keys, valid_array) if keep]
    return selected_keys, labels.loc[valid].tolist(), invalid_labels


def _deduplicate_records(
    keys: Sequence[Any],
    values: Sequence[Any],
    source: str,
    conflict_policy: str,
) -> tuple[dict[Any, Any], int, int]:
    mapping: dict[Any, Any] = {}
    duplicate_count = 0
    conflict_count = 0
    for key, value in zip(keys, values):
        if key not in mapping:
            mapping[key] = value
            continue
        duplicate_count += 1
        if str(mapping[key]) == str(value):
            continue
        conflict_count += 1
        if conflict_policy == "error":
            raise ValueError(
                f"Conflicting annotations for key {key!r} in {source}: "
                f"{mapping[key]!r} versus {value!r}."
            )
        if conflict_policy == "last":
            mapping[key] = value
    return mapping, duplicate_count, conflict_count


def overlay_annotations(
    input_paths: Sequence[str],
    *,
    table_key: str = "",
    annotation_csv: str = "",
    annotation_col: str = "auto",
    target_col: str = "sub_celltype",
    fallback_col: str = "celltype",
    csv_cell_col: str = "auto",
    csv_region_col: str = "auto",
    input_cell_col: str = "auto",
    input_region_col: str = "auto",
    conflict_policy: str = "error",
    existing_policy: str = "overwrite_matched",
    min_match_rate: float = 0.95,
) -> tuple[sd.SpatialData, str, list[dict[str, Any]]]:
    if not input_paths:
        raise ValueError("merge_by=annotation requires a parent SpatialData input.")
    csv_paths = _parse_csv_inputs(annotation_csv)
    zarr_paths = list(input_paths[1:])
    if csv_paths and zarr_paths:
        raise ValueError("Annotation merge cannot mix CSV and Zarr sources in one run.")
    if not csv_paths and not zarr_paths:
        raise ValueError(
            "Annotation merge requires --annotation_csv or one or more subset Zarr inputs after the parent."
        )
    if csv_paths and len(input_paths) != 1:
        raise ValueError("CSV annotation mode accepts exactly one INPUT: the parent Zarr.")

    parent_path = input_paths[0]
    parent = _read_sdata(parent_path)
    parent_table_key = _resolve_table_key(parent, table_key, parent_path)
    _validate_table_links(parent, parent_table_key, parent_path)
    original = parent.tables[parent_table_key]
    table = original.copy()
    context = _base_match_context(
        table,
        parent_path,
        input_cell_col,
        input_region_col,
    )
    base_key_set = set(context["keys"])

    sources = csv_paths if csv_paths else zarr_paths
    source_kind = "csv" if csv_paths else "zarr"
    combined: dict[Any, Any] = {}
    report: list[dict[str, Any]] = []

    for source_path in sources:
        if source_kind == "csv":
            keys, values, invalid_labels = _records_from_csv(
                source_path,
                annotation_col=annotation_col,
                csv_cell_col=csv_cell_col,
                csv_region_col=csv_region_col,
                use_region=context["use_region"],
            )
        else:
            keys, values, invalid_labels = _records_from_zarr(
                source_path,
                table_key=table_key,
                annotation_col=annotation_col,
                use_region=context["use_region"],
            )
        local, duplicate_count, conflict_count = _deduplicate_records(
            keys, values, source_path, conflict_policy
        )
        requested = len(local)
        matched_keys = set(local).intersection(base_key_set)
        matched = len(matched_keys)
        if requested == 0 or matched == 0:
            raise ValueError(f"No valid annotation IDs from {source_path} matched the parent table.")
        match_rate = matched / requested
        if match_rate < min_match_rate:
            raise ValueError(
                f"Annotation match rate for {source_path} is {match_rate:.3f}, below "
                f"min_match_rate={min_match_rate:.3f}."
            )

        cross_conflicts = 0
        for match_key in matched_keys:
            value = local[match_key]
            if match_key not in combined:
                combined[match_key] = value
                continue
            if str(combined[match_key]) == str(value):
                continue
            cross_conflicts += 1
            if conflict_policy == "error":
                raise ValueError(
                    f"Conflicting annotations across sources for key {match_key!r}: "
                    f"{combined[match_key]!r} versus {value!r}."
                )
            if conflict_policy == "last":
                combined[match_key] = value

        report.append(
            {
                "mode": "annotation",
                "source": source_path,
                "source_type": source_kind,
                "requested": requested,
                "matched": matched,
                "unmatched": requested - matched,
                "duplicate_ids": duplicate_count,
                "conflicts": conflict_count + cross_conflicts,
                "invalid_labels": invalid_labels,
                "match_rate": match_rate,
                "target_col": target_col,
            }
        )

    target_existed = target_col in table.obs
    if target_existed and existing_policy == "error":
        raise ValueError(
            f"Target column {target_col!r} already exists and existing_policy='error'."
        )
    if target_existed:
        merged_labels = table.obs[target_col].astype(object).copy()
    elif not _is_auto(fallback_col):
        if fallback_col not in table.obs:
            raise ValueError(f"Fallback column {fallback_col!r} was not found in the parent table.")
        merged_labels = table.obs[fallback_col].astype(object).copy()
    else:
        merged_labels = pd.Series(pd.NA, index=table.obs_names, dtype=object)

    key_to_position = {key: position for position, key in enumerate(context["keys"])}
    updated = 0
    for match_key, value in combined.items():
        position = key_to_position[match_key]
        if target_existed and existing_policy == "fill_missing":
            current = merged_labels.iloc[position]
            if pd.notna(current) and str(current).strip() != "":
                continue
        merged_labels.iloc[position] = value
        updated += 1
    if updated == 0:
        raise ValueError("No parent annotations were updated under the selected existing_policy.")

    table.obs[target_col] = pd.Categorical(merged_labels)
    _append_merge_history(
        table.uns,
        {
            "mode": "annotation",
            "sources": list(sources),
            "target_col": target_col,
            "fallback_col": None if _is_auto(fallback_col) else fallback_col,
            "updated": updated,
            "conflict_policy": conflict_policy,
            "existing_policy": existing_policy,
        },
    )

    if table.n_obs != original.n_obs or table.n_vars != original.n_vars:
        raise RuntimeError("Annotation overlay changed the parent table dimensions.")
    if not table.obs_names.equals(original.obs_names):
        raise RuntimeError("Annotation overlay changed parent observation IDs or order.")
    if not table.var_names.equals(original.var_names):
        raise RuntimeError("Annotation overlay changed parent feature IDs or order.")
    if set(table.layers) != set(original.layers) or set(table.obsm) != set(original.obsm):
        raise RuntimeError("Annotation overlay changed parent layers or obsm keys.")

    parent.tables[parent_table_key] = table
    _validate_table_links(parent, parent_table_key, "annotation output")
    return parent, parent_table_key, report


def _write_validated(
    sdata: sd.SpatialData,
    output_path: str,
    table_key: str,
    expected_obs: int,
) -> None:
    output_path = os.path.abspath(output_path)
    output_parent = os.path.dirname(output_path)
    os.makedirs(output_parent, exist_ok=True)
    temporary = os.path.join(
        output_parent,
        f".{os.path.basename(output_path)}.tmp-{uuid.uuid4().hex}",
    )
    backup = os.path.join(
        output_parent,
        f".{os.path.basename(output_path)}.backup-{uuid.uuid4().hex}",
    )
    try:
        sdata.write(temporary, overwrite=False)
        reopened = sd.read_zarr(temporary)
        if table_key not in reopened.tables:
            raise RuntimeError(f"Reopened output is missing table {table_key!r}.")
        if reopened.tables[table_key].n_obs != expected_obs:
            raise RuntimeError("Reopened output has an unexpected observation count.")
        _validate_table_links(reopened, table_key, "reopened output")

        had_existing = os.path.exists(output_path)
        if had_existing:
            os.replace(output_path, backup)
        try:
            os.replace(temporary, output_path)
        except Exception:
            if had_existing and os.path.exists(backup):
                os.replace(backup, output_path)
            raise
        if had_existing and os.path.exists(backup):
            shutil.rmtree(backup)
    finally:
        if os.path.exists(temporary):
            shutil.rmtree(temporary)
        if os.path.exists(backup) and os.path.exists(output_path):
            shutil.rmtree(backup)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge SpatialData samples or annotations")
    parser.add_argument("--INPUT", nargs="+", required=True, help="Input SpatialData Zarr directories")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--output_name", default="concatenated_sdata.zarr", help="Output Zarr name")
    parser.add_argument("--merge_by", required=True, type=str.lower, choices=("sample", "annotation"))
    parser.add_argument("--table_key", default="", help="Table to merge; empty requires an unambiguous table")

    parser.add_argument("--sample_ids", default="", help="Comma-separated source IDs in input order")
    parser.add_argument("--sample_col", default="sample", help="Observation sample column")
    parser.add_argument("--feature_join", default="inner", type=str.lower, choices=("inner", "outer"))
    parser.add_argument("--label_cols", default="", help="Comma-separated label columns")
    parser.add_argument(
        "--label_policy",
        default="preserve",
        type=str.lower,
        choices=("preserve", "prefix", "offset"),
    )

    parser.add_argument("--annotation_csv", default="", help="CSV file, directory, or comma-separated CSV paths")
    parser.add_argument("--annotation_col", default="auto", help="Annotation column in CSV or subset Zarr")
    parser.add_argument("--target_col", default="sub_celltype", help="Parent obs column to update")
    parser.add_argument("--fallback_col", default="celltype", help="Parent column used to initialize target_col")
    parser.add_argument("--csv_cell_col", default="auto", help="Cell ID column in annotation CSV")
    parser.add_argument("--csv_region_col", default="auto", help="Region/sample column in annotation CSV")
    parser.add_argument("--input_cell_col", default="auto", help="Cell ID column in parent table")
    parser.add_argument("--input_region_col", default="auto", help="Region column in parent table")
    parser.add_argument(
        "--conflict_policy",
        default="error",
        type=str.lower,
        choices=("error", "first", "last"),
    )
    parser.add_argument(
        "--existing_policy",
        default="overwrite_matched",
        type=str.lower,
        choices=("overwrite_matched", "fill_missing", "error"),
    )
    parser.add_argument("--min_match_rate", type=float, default=0.95)
    return parser


def main(argv: Sequence[str] | None = None) -> str:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not 0.0 <= args.min_match_rate <= 1.0:
        parser.error("--min_match_rate must be between 0 and 1")

    log_step(logger, 1, 3, f"loading {len(args.INPUT)} input file(s)")
    if args.merge_by == "sample":
        result, table_key, report = merge_samples(
            args.INPUT,
            table_key=args.table_key,
            sample_ids=args.sample_ids,
            sample_col=args.sample_col,
            feature_join=args.feature_join,
            label_cols=args.label_cols,
            label_policy=args.label_policy,
        )
    else:
        result, table_key, report = overlay_annotations(
            args.INPUT,
            table_key=args.table_key,
            annotation_csv=args.annotation_csv,
            annotation_col=args.annotation_col,
            target_col=args.target_col,
            fallback_col=args.fallback_col,
            csv_cell_col=args.csv_cell_col,
            csv_region_col=args.csv_region_col,
            input_cell_col=args.input_cell_col,
            input_region_col=args.input_region_col,
            conflict_policy=args.conflict_policy,
            existing_policy=args.existing_policy,
            min_match_rate=args.min_match_rate,
        )

    log_step(logger, 2, 3, f"writing {args.merge_by} merge output")
    output_path = os.path.join(args.output_dir, args.output_name)
    expected_obs = result.tables[table_key].n_obs
    _write_validated(result, output_path, table_key, expected_obs)
    pd.DataFrame(report).to_csv(os.path.join(args.output_dir, "merge_report.csv"), index=False)
    log_step(logger, 3, 3, f"merged output saved to {output_path}")
    return output_path


if __name__ == "__main__":
    main()
