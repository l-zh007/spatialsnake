import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as spd
from scipy import sparse

from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("RCTD_validate_input")


def _validated_names(index, axis_name):
    """Return string identifiers after rejecting missing, blank, or duplicate names."""

    names = pd.Index(index)
    if names.hasnans:
        raise ValueError(f"{axis_name} contains missing identifiers.")

    names = pd.Index(names.astype(str))
    blank = names.str.strip() == ""
    if blank.any():
        raise ValueError(f"{axis_name} contains {int(blank.sum())} blank identifier(s).")
    if not names.is_unique:
        duplicates = names[names.duplicated(keep=False)].unique().tolist()[:10]
        raise ValueError(f"{axis_name} contains duplicate identifiers: {duplicates}")
    return names


def _validate_adata_axes(adata, object_name):
    if adata.n_obs == 0:
        raise ValueError(f"{object_name} contains no observations.")
    if adata.n_vars == 0:
        raise ValueError(f"{object_name} contains no genes.")

    obs_names = _validated_names(adata.obs_names, f"{object_name}.obs_names")
    var_names = _validated_names(adata.var_names, f"{object_name}.var_names")
    return obs_names, var_names


def _validate_spatial_coordinates(adata):
    if "spatial" not in adata.obsm:
        raise ValueError("Converted H5AD does not contain obsm['spatial'] coordinates.")

    coordinates = np.asarray(adata.obsm["spatial"])
    expected = f"({adata.n_obs}, >=2)"
    if (
        coordinates.ndim != 2
        or coordinates.shape[0] != adata.n_obs
        or coordinates.shape[1] < 2
    ):
        raise ValueError(
            f"Converted H5AD has invalid obsm['spatial'] shape {coordinates.shape}; "
            f"expected {expected}."
        )

    try:
        coordinates = np.asarray(coordinates[:, :2], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Converted H5AD spatial coordinates must be numeric.") from exc
    if not np.isfinite(coordinates).all():
        raise ValueError("Converted H5AD spatial coordinates contain non-finite values.")


def _count_matrix_error(matrix, expected_shape):
    if matrix is None:
        return "is missing"
    if tuple(matrix.shape) != tuple(expected_shape):
        return f"shape {matrix.shape} does not match {expected_shape}"

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size == 0:
        return "contains no positive counts"

    try:
        finite = np.isfinite(values)
    except TypeError:
        return "contains non-numeric values"
    if not finite.all():
        return "contains non-finite values"
    if np.any(values < 0):
        return "contains negative values"
    if not np.equal(values, np.floor(values)).all():
        return "contains non-integer values"
    if not np.any(values > 0):
        return "contains no positive counts"
    return None


def _validate_count_source(adata, var_names):
    candidates = [("X", adata.X)]
    for layer_name in ("counts", "raw_counts"):
        if layer_name in adata.layers:
            candidates.append((f"layers['{layer_name}']", adata.layers[layer_name]))

    candidate_errors = []
    raw_is_aligned = False
    if adata.raw is not None:
        raw_var_names = _validated_names(
            adata.raw.var_names, "Converted H5AD.raw.var_names"
        )
        if adata.raw.n_obs != adata.n_obs:
            candidate_errors.append(
                f"raw.X: observation count {adata.raw.n_obs} does not match {adata.n_obs}"
            )
        elif not var_names.isin(raw_var_names).all():
            missing = var_names[~var_names.isin(raw_var_names)].tolist()[:10]
            candidate_errors.append(
                "raw.X: genes cannot be aligned to converted var_names; "
                f"missing examples={missing}"
            )
        else:
            raw_is_aligned = True

    for source_name, matrix in candidates:
        error = _count_matrix_error(matrix, adata.shape)
        if error is None:
            logger.info("Using validated non-negative integer counts from %s", source_name)
            return source_name, matrix, adata.var.copy()
        candidate_errors.append(f"{source_name}: {error}")

    if raw_is_aligned:
        raw_error = _count_matrix_error(
            adata.raw.X, (adata.raw.n_obs, adata.raw.n_vars)
        )
        if raw_error is None:
            logger.info(
                "Using validated non-negative integer counts from raw.X "
                "(gene-aligned)"
            )
            if sparse.issparse(adata.raw.X):
                return "raw.X (gene-aligned)", None, None
            return (
                "raw.X (gene-aligned)",
                adata.raw.X,
                adata.raw.var.copy(),
            )
        candidate_errors.append(f"raw.X (gene-aligned): {raw_error}")

    details = "; ".join(candidate_errors)
    raise ValueError(
        "Converted H5AD does not contain usable non-negative integer counts in X, "
        "layers['counts'], layers['raw_counts'], or gene-aligned raw.X. "
        f"Checked candidates: {details}"
    )


def validate_rctd_input(spatial_zarr, converted_h5ad):
    spatial_zarr = os.path.normpath(os.path.abspath(spatial_zarr))
    converted_h5ad = os.path.abspath(converted_h5ad)

    if not spatial_zarr.lower().endswith(".zarr"):
        raise ValueError(
            f"RCTD spatial input must be a .zarr directory: {spatial_zarr}"
        )
    if not os.path.isdir(spatial_zarr):
        raise FileNotFoundError(
            f"RCTD SpatialData Zarr directory does not exist: {spatial_zarr}"
        )
    if not os.path.isfile(converted_h5ad):
        raise FileNotFoundError(
            f"Converted RCTD H5AD file does not exist: {converted_h5ad}"
        )
    if os.path.getsize(converted_h5ad) == 0:
        raise ValueError(f"Converted RCTD H5AD file is empty: {converted_h5ad}")

    log_step(logger, 1, 4, f"reading source SpatialData: {spatial_zarr}")
    source_sdata = spd.read_zarr(spatial_zarr)
    table_keys = list(source_sdata.tables.keys())
    if not table_keys:
        raise ValueError("Source SpatialData Zarr does not contain any table.")
    if len(table_keys) != 1:
        raise ValueError(
            "RCTD requires exactly one SpatialData table so conversion, fitting, "
            f"and result write-back use the same observations; found: {table_keys}"
        )
    table_key = table_keys[0]
    source_table = source_sdata.tables[table_key]
    source_obs_names, _ = _validate_adata_axes(
        source_table, f"Source SpatialData table '{table_key}'"
    )

    log_step(logger, 2, 4, f"reading converted H5AD: {converted_h5ad}")
    try:
        converted = ad.read_h5ad(converted_h5ad)
    except Exception as exc:
        raise ValueError(
            f"Converted RCTD H5AD could not be read: {converted_h5ad}"
        ) from exc
    converted_obs_names, converted_var_names = _validate_adata_axes(
        converted, "Converted H5AD"
    )

    source_ids = set(source_obs_names)
    converted_ids = set(converted_obs_names)
    if source_ids != converted_ids:
        missing = sorted(source_ids.difference(converted_ids))[:10]
        extra = sorted(converted_ids.difference(source_ids))[:10]
        raise ValueError(
            "Converted H5AD observation identifiers do not match the source SpatialData "
            f"table '{table_key}'. Missing from H5AD ({len(source_ids - converted_ids)}): "
            f"{missing}; extra in H5AD ({len(converted_ids - source_ids)}): {extra}"
        )
    if not source_obs_names.equals(converted_obs_names):
        logger.warning(
            "Converted H5AD contains the same observation identifiers in a different order"
        )

    log_step(logger, 3, 4, "validating spatial coordinates and count matrices")
    _validate_spatial_coordinates(converted)
    count_source, count_matrix, count_var = _validate_count_source(
        converted, converted_var_names
    )
    if count_matrix is not None:
        # RCTD.R deliberately imports the spatial H5AD with schard's
        # ``use.raw=TRUE``.  Materialize the validated source in raw/X so that
        # the matrix checked here is exactly the one used for model fitting.
        raw_counts = (
            count_matrix.tocsr(copy=True)
            if sparse.issparse(count_matrix)
            else sparse.csr_matrix(np.asarray(count_matrix))
        )
        counts_adata = ad.AnnData(
            X=raw_counts,
            obs=pd.DataFrame(index=converted.obs_names.copy()),
            var=count_var,
        )
        converted.raw = counts_adata
        converted.write_h5ad(converted_h5ad)
        logger.info(
            "Materialized validated %s counts in temporary H5AD raw.X",
            count_source,
        )

    log_step(
        logger,
        4,
        4,
        (
            f"validated {converted.n_obs} observations and {converted.n_vars} genes "
            f"(counts source: {count_source})"
        ),
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate the temporary H5AD generated for an RCTD run."
    )
    parser.add_argument(
        "--spatial_zarr",
        required=True,
        help="Original SpatialData .zarr directory from sample.txt.",
    )
    parser.add_argument(
        "--converted_h5ad",
        required=True,
        help="Temporary H5AD generated from the SpatialData Zarr.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    validate_rctd_input(args.spatial_zarr, args.converted_h5ad)
