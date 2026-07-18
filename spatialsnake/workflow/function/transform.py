"""Lightweight SpatialData/AnnData/Seurat format conversion.

The public entry point remains ``spatialsnake useful_tool --option=transform``.
SpatialData is exported through the project's point-compatible legacy AnnData
converter.  Schard is used only for the one-way H5AD-to-Seurat step.
"""

from __future__ import annotations

import argparse
import copy
import gc
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Iterable, Sequence

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import spatialdata as sd
from spatialdata.models import TableModel
from spatialdata_io.experimental import from_legacy_anndata

from spatialsnake.workflow.function.legacy_anndata import to_legacy_anndata
from spatialsnake.workflow.function.logging_utils import log_step, setup_logger


logger = setup_logger("useful_transform")
GIB = 1024**3
SCHARD_MIN_VERSION = "1.1.0"


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y", "t"}:
        return True
    if normalized in {"false", "0", "no", "n", "f", "", "none", "null"}:
        return False
    raise ValueError(f"Expected a boolean value, got {value!r}")


def _clean_stem(path: str | Path) -> str:
    name = Path(path).name.rstrip("/")
    return name[:-5] if name.lower().endswith(".zarr") else Path(name).stem


def _select_table(sdata: sd.SpatialData, table_key: str = ""):
    if table_key:
        if table_key not in sdata.tables:
            raise ValueError(
                f"table_key={table_key!r} is not present; available tables: {list(sdata.tables)}"
            )
        return table_key, sdata.tables[table_key]
    if len(sdata.tables) != 1:
        raise ValueError(
            "Transform requires exactly one unambiguous table. "
            f"Found {list(sdata.tables)}; set --table_key."
        )
    key = next(iter(sdata.tables))
    return key, sdata.tables[key]


def _spatial_attrs(table) -> dict[str, Any]:
    attrs = table.uns.get(TableModel.ATTRS_KEY, {})
    return attrs if isinstance(attrs, dict) else {}


def _table_regions(table) -> list[str]:
    region = _spatial_attrs(table).get("region")
    if region is None:
        return []
    values = region if isinstance(region, (list, tuple, np.ndarray, pd.Index)) else [region]
    return [str(value) for value in values]


def _coordinate_system_for_region(sdata: sd.SpatialData, region: str) -> str:
    systems = sorted(str(value) for value in sdata.coordinate_systems)
    if region in systems:
        return region

    candidates: list[str] = []
    for system in systems:
        filtered = sdata.filter_by_coordinate_system(system)
        if region in filtered.shapes or region in filtered.labels or region in filtered.points:
            candidates.append(system)
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError(
            f"Region {region!r} is not connected to any coordinate system; available: {systems}"
        )
    raise ValueError(
        f"Region {region!r} is connected to multiple coordinate systems {candidates}; "
        "the conversion is ambiguous."
    )


def _unique_library_id(candidate: str, source: str, used: set[str]) -> str:
    value = candidate
    if value in used:
        value = f"{candidate}-{source}"
    counter = 2
    base = value
    while value in used:
        value = f"{base}-{counter}"
        counter += 1
    used.add(value)
    return value


def _normalize_spatial_uns(adata: ad.AnnData, library_id: str) -> list[str]:
    warnings: list[str] = []
    spatial = adata.uns.get("spatial")
    if not isinstance(spatial, dict) or not spatial:
        return warnings

    if library_id in spatial:
        entry = spatial[library_id]
    elif len(spatial) == 1:
        entry = next(iter(spatial.values()))
    else:
        keys = sorted(spatial)
        image_keys = [key for key in keys if isinstance(spatial[key], dict) and spatial[key].get("images")]
        chosen = image_keys[0] if image_keys else keys[0]
        entry = spatial[chosen]
        warnings.append(
            f"multiple legacy image entries {keys} were available for {library_id}; used {chosen}"
        )
    adata.uns["spatial"] = {library_id: entry}
    return warnings


def _count_spatial_images(adata: ad.AnnData) -> int:
    spatial = adata.uns.get("spatial")
    if not isinstance(spatial, dict):
        return 0
    return sum(
        len(entry.get("images", {}))
        for entry in spatial.values()
        if isinstance(entry, dict) and isinstance(entry.get("images", {}), dict)
    )


def _validate_spatial_coordinates(adata: ad.AnnData, context: str) -> None:
    if "spatial" not in adata.obsm:
        return
    coords = np.asarray(adata.obsm["spatial"])
    if coords.ndim != 2 or coords.shape[0] != adata.n_obs or coords.shape[1] < 2:
        raise ValueError(f"{context}: obsm['spatial'] must have shape (n_obs, >=2); got {coords.shape}")
    if not np.isfinite(coords[:, :2]).all():
        raise ValueError(f"{context}: obsm['spatial'] contains missing or non-finite coordinates")


def _export_library(
    sdata: sd.SpatialData,
    table_key: str,
    region: str | None,
    library_id: str,
    save_image: bool,
) -> tuple[ad.AnnData, dict[str, Any]]:
    if region is None:
        adata = sdata.tables[table_key].copy()
        coordinate_system = ""
    else:
        coordinate_system = _coordinate_system_for_region(sdata, region)
        adata = to_legacy_anndata(
            sdata,
            table_name=table_key,
            coordinate_system=coordinate_system,
            include_images=save_image,
        )

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"Library {library_id!r} produced an empty AnnData object")
    if not adata.obs_names.is_unique:
        raise ValueError(f"Library {library_id!r} contains duplicated observation IDs")
    _validate_spatial_coordinates(adata, library_id)
    messages = _normalize_spatial_uns(adata, library_id)
    adata.obs["library_id"] = pd.Categorical([library_id] * adata.n_obs)
    image_count = _count_spatial_images(adata) if save_image else 0
    return adata, {
        "library_id": library_id,
        "region": region or "",
        "coordinate_system": coordinate_system,
        "observations": adata.n_obs,
        "genes": adata.n_vars,
        "images": image_count,
        "warnings": "; ".join(messages),
    }


def _resolve_duplicate_obs_names(paths: Sequence[Path], library_ids: Sequence[str]) -> bool:
    all_names: list[str] = []
    per_file: list[pd.Index] = []
    for path in paths:
        backed = ad.read_h5ad(path, backed="r")
        names = pd.Index(backed.obs_names.astype(str))
        backed.file.close()
        if not names.is_unique:
            raise ValueError(f"Duplicated observation IDs within {path}")
        per_file.append(names)
        all_names.extend(names.tolist())

    duplicated = pd.Index(all_names).duplicated(keep=False)
    duplicate_values = set(pd.Index(all_names)[duplicated])
    if not duplicate_values:
        return False

    for path, library_id, names in zip(paths, library_ids, per_file, strict=True):
        renamed = [
            name if name not in duplicate_values or name.endswith(f"-{library_id}") else f"{name}-{library_id}"
            for name in names
        ]
        if renamed == names.tolist():
            continue
        adata = ad.read_h5ad(path)
        adata.obs_names = pd.Index(renamed, name=adata.obs_names.name)
        for key in adata.obsm:
            if isinstance(adata.obsm[key], pd.DataFrame):
                frame = adata.obsm[key].copy()
                frame.index = adata.obs_names
                adata.obsm[key] = frame
        adata.write_h5ad(path)

    final_names: list[str] = []
    for path in paths:
        backed = ad.read_h5ad(path, backed="r")
        final_names.extend(backed.obs_names.astype(str).tolist())
        backed.file.close()
    if not pd.Index(final_names).is_unique:
        raise ValueError("Observation IDs remain duplicated after adding library suffixes")
    return True


def _replace_file(staged: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        shutil.rmtree(destination) if destination.is_dir() else destination.unlink()
    os.replace(staged, destination)


def _restore_pairwise_matrices(paths: Sequence[Path], destination: Path) -> None:
    """Append block-diagonal ``obsp`` matrices after on-disk concatenation.

    AnnData 0.12 does not yet implement ``pairwise=True`` for
    ``concat_on_disk``.  Pairwise graphs are therefore assembled one key at a
    time as sparse block diagonals and written directly to the output store.
    This preserves within-library graphs without introducing cross-library
    edges or materialising the expression matrix.
    """
    backed_inputs = [ad.read_h5ad(path, backed="r") for path in paths]
    try:
        keys = sorted({key for backed in backed_inputs for key in backed.obsp})
        if not keys:
            return
        with h5py.File(destination, "r+") as handle:
            group = handle.require_group("obsp")
            group.attrs["encoding-type"] = "dict"
            group.attrs["encoding-version"] = "0.1.0"
            for key in keys:
                blocks = []
                for backed in backed_inputs:
                    if key not in backed.obsp:
                        blocks.append(sparse.csr_matrix((backed.n_obs, backed.n_obs)))
                        continue
                    value = backed.obsp[key]
                    value = value.to_memory() if hasattr(value, "to_memory") else value
                    blocks.append(value if sparse.issparse(value) else sparse.csr_matrix(value))
                matrix = sparse.block_diag(blocks, format="csr")
                if key in group:
                    del group[key]
                ad.io.write_elem(group, key, matrix)
    finally:
        for backed in backed_inputs:
            backed.file.close()


def _restore_var_metadata(paths: Sequence[Path], destination: Path) -> None:
    """Merge variable metadata by gene name after matrix-only concatenation.

    AnnData 0.12 can produce non-serialisable pandas ``Series`` values when
    ``concat_on_disk(..., merge='first')`` handles real-world ``var`` columns.
    The matrix join remains on disk; only the comparatively small gene metadata
    frames are loaded and combined here.
    """
    output = ad.read_h5ad(destination, backed="r")
    target_index = pd.Index(output.var_names.astype(str), name=output.var_names.name)
    output.file.close()

    frames: list[pd.DataFrame] = []
    all_columns: list[str] = []
    for path in paths:
        backed = ad.read_h5ad(path, backed="r")
        frame = backed.var.copy()
        backed.file.close()
        frame.index = frame.index.astype(str)
        frames.append(frame)
        all_columns.extend(str(column) for column in frame.columns)

    merged = pd.DataFrame(index=target_index)
    for column in dict.fromkeys(all_columns):
        combined = pd.Series(index=target_index, dtype=object)
        conflicts = 0
        for frame in frames:
            if column not in frame:
                continue
            values = frame[column]
            if isinstance(values.dtype, pd.CategoricalDtype):
                values = values.astype(object)
            values = values.reindex(target_index)
            overlap = combined.notna() & values.notna()
            if overlap.any():
                conflicts += int((combined.loc[overlap].astype(str) != values.loc[overlap].astype(str)).sum())
            fill = combined.isna() & values.notna()
            if fill.any():
                combined.loc[fill] = values.loc[fill].astype(object)
        if conflicts:
            logger.warning(
                "Variable metadata column %s contained %d conflicting values; kept the first value",
                column,
                conflicts,
            )
        merged[column] = combined.infer_objects(copy=False)

    with h5py.File(destination, "r+") as handle:
        if "var" in handle:
            del handle["var"]
        ad.io.write_elem(handle, "var", merged)


def _write_raw_h5ad_view(source: Path, destination: Path) -> None:
    """Create a lightweight H5AD view whose X/var point to the source raw group."""
    with h5py.File(source, "r") as source_handle, h5py.File(destination, "w") as target:
        if "raw" not in source_handle:
            raise ValueError(f"Raw matrix is missing from {source}")
        for key, value in source_handle.attrs.items():
            target.attrs[key] = value
        source_name = str(source.resolve())
        target["X"] = h5py.ExternalLink(source_name, "/raw/X")
        target["obs"] = h5py.ExternalLink(source_name, "/obs")
        target["var"] = h5py.ExternalLink(source_name, "/raw/var")
        if "/raw/varm" in source_handle:
            target["varm"] = h5py.ExternalLink(source_name, "/raw/varm")
        for key in ("obsm", "obsp", "varm", "varp", "layers", "uns"):
            if key in target:
                continue
            group = target.create_group(key)
            group.attrs["encoding-type"] = "dict"
            group.attrs["encoding-version"] = "0.1.0"


def _restore_raw_matrix(paths: Sequence[Path], destination: Path) -> None:
    raw_presence = []
    for path in paths:
        with h5py.File(path, "r") as handle:
            raw_presence.append("raw" in handle)
    if not any(raw_presence):
        return
    if not all(raw_presence):
        raise ValueError(
            "Cannot concatenate inputs when only some libraries contain .raw. "
            "Provide consistent Spatialsnake tables or convert libraries separately."
        )

    with tempfile.TemporaryDirectory(prefix="spatialsnake-raw-", dir=destination.parent) as temp_dir:
        raw_views = []
        for index, path in enumerate(paths):
            view = Path(temp_dir) / f"raw-{index}.h5ad"
            _write_raw_h5ad_view(path, view)
            raw_views.append(view)
        combined = Path(temp_dir) / "raw-combined.h5ad"
        ad.experimental.concat_on_disk(
            [str(path) for path in raw_views],
            str(combined),
            join="outer",
            merge=None,
            uns_merge="first",
            fill_value=0,
            max_loaded_elems=25_000_000,
        )
        _restore_var_metadata(raw_views, combined)
        with h5py.File(paths[0], "r") as first, h5py.File(combined, "r") as raw_source, h5py.File(
            destination, "r+"
        ) as target:
            if "raw" in target:
                del target["raw"]
            raw_group = target.create_group("raw")
            for key, value in first["raw"].attrs.items():
                raw_group.attrs[key] = value
            for key in ("X", "var", "varm"):
                if key in raw_source:
                    raw_source.copy(key, raw_group)


def _restore_spatial_uns(paths: Sequence[Path], destination: Path) -> None:
    """Merge per-library legacy image metadata without loading image arrays."""
    sources_with_spatial = []
    for path in paths:
        with h5py.File(path, "r") as handle:
            spatial = handle.get("/uns/spatial")
            if isinstance(spatial, h5py.Group) and len(spatial):
                sources_with_spatial.append(path)
    if not sources_with_spatial:
        return

    with h5py.File(destination, "r+") as target:
        uns = target.require_group("uns")
        uns.attrs["encoding-type"] = "dict"
        uns.attrs["encoding-version"] = "0.1.0"
        if "spatial" in uns:
            del uns["spatial"]
        merged = uns.create_group("spatial")
        merged.attrs["encoding-type"] = "dict"
        merged.attrs["encoding-version"] = "0.1.0"
        for path in sources_with_spatial:
            with h5py.File(path, "r") as source:
                source_spatial = source["/uns/spatial"]
                for library_id in source_spatial:
                    if library_id in merged:
                        raise ValueError(
                            f"Duplicate uns['spatial'] library ID {library_id!r} while concatenating"
                        )
                    source.copy(source_spatial[library_id], merged, name=library_id)


def _concat_h5ad_on_disk(paths: Sequence[Path], destination: Path) -> None:
    staged = destination.parent / f".{destination.stem}.tmp-{os.getpid()}.h5ad"
    if staged.exists():
        staged.unlink()
    if len(paths) == 1:
        shutil.copy2(paths[0], staged)
    else:
        ad.experimental.concat_on_disk(
            [str(path) for path in paths],
            str(staged),
            join="outer",
            merge=None,
            uns_merge="unique",
            fill_value=0,
            max_loaded_elems=25_000_000,
        )
        _restore_var_metadata(paths, staged)
        _restore_raw_matrix(paths, staged)
        _restore_spatial_uns(paths, staged)
        _restore_pairwise_matrices(paths, staged)
    _replace_file(staged, destination)


def _validate_h5ad(path: str | Path, expected_obs: int | None = None) -> dict[str, Any]:
    backed = ad.read_h5ad(path, backed="r")
    try:
        if expected_obs is not None and backed.n_obs != expected_obs:
            raise ValueError(f"Expected {expected_obs} observations, found {backed.n_obs}")
        if backed.n_obs == 0 or backed.n_vars == 0:
            raise ValueError("Converted H5AD is empty")
        if not backed.obs_names.is_unique:
            raise ValueError("Converted H5AD contains duplicated observation IDs")
        if backed.raw is not None and backed.raw.n_obs != backed.n_obs:
            raise ValueError("Converted H5AD raw.X is not aligned to observations")
        if "spatial" in backed.obsm:
            _validate_spatial_coordinates(backed, str(path))
        libraries = (
            sorted(backed.obs["library_id"].astype(str).unique())
            if "library_id" in backed.obs
            else []
        )
        return {
            "observations": backed.n_obs,
            "genes": backed.n_vars,
            "libraries": libraries,
            "has_raw": backed.raw is not None,
        }
    finally:
        backed.file.close()


def _export_zarr_to_h5ad(
    input_list: Sequence[str],
    output_h5ad: str | Path,
    save_image: bool = True,
    table_key: str = "",
) -> pd.DataFrame:
    if not input_list:
        raise ValueError("At least one Zarr input is required")
    output = Path(output_h5ad)
    output.parent.mkdir(parents=True, exist_ok=True)
    used_library_ids: set[str] = set()
    records: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="spatialsnake-transform-", dir=output.parent) as temp_dir:
        temp_paths: list[Path] = []
        library_ids: list[str] = []
        for input_index, file_path in enumerate(input_list):
            logger.info("Reading SpatialData Zarr: %s", file_path)
            sdata = sd.read_zarr(file_path)
            selected_key, table = _select_table(sdata, table_key)
            regions = _table_regions(table) or [None]
            source = _clean_stem(file_path)
            for region_index, region in enumerate(regions):
                proposed = region or source
                library_id = _unique_library_id(proposed, source, used_library_ids)
                logger.info(
                    "Exporting input=%s region=%s as library_id=%s",
                    file_path,
                    region or "table-only",
                    library_id,
                )
                library, record = _export_library(
                    sdata,
                    selected_key,
                    region,
                    library_id,
                    save_image,
                )
                record.update({"input": str(file_path), "table_key": selected_key})
                path = Path(temp_dir) / f"library-{input_index}-{region_index}.h5ad"
                library.write_h5ad(path)
                temp_paths.append(path)
                library_ids.append(library_id)
                records.append(record)
                del library
                gc.collect()

        renamed = _resolve_duplicate_obs_names(temp_paths, library_ids)
        _concat_h5ad_on_disk(temp_paths, output)

    expected_obs = sum(int(record["observations"]) for record in records)
    summary = _validate_h5ad(output, expected_obs)
    logger.info(
        "Validated H5AD: obs=%d genes=%d libraries=%s; duplicate IDs renamed=%s",
        summary["observations"],
        summary["genes"],
        summary["libraries"],
        renamed,
    )
    return pd.DataFrame(records)


def zarr_to_h5ad(
    INPUT_list,
    output_h5ad,
    save_image=True,
    table_key="",
    return_adata=True,
):
    """Backward-compatible public wrapper used by existing workflow scripts."""
    _export_zarr_to_h5ad(
        [str(value) for value in INPUT_list],
        output_h5ad,
        save_image=parse_bool(save_image),
        table_key=table_key,
    )
    return ad.read_h5ad(output_h5ad) if return_adata else None


def _h5_group_bytes(group: h5py.Group | h5py.Dataset) -> int:
    if isinstance(group, h5py.Dataset):
        return int(np.prod(group.shape, dtype=np.int64)) * int(group.dtype.itemsize)
    return sum(_h5_group_bytes(value) for value in group.values())


def _h5_matrix_bytes(path: str | Path, source: str) -> int:
    key = "/raw/X" if source == "raw" else "/X"
    with h5py.File(path, "r") as handle:
        if key not in handle:
            return 0
        return _h5_group_bytes(handle[key])


def _h5_spatial_bytes(path: str | Path) -> int:
    with h5py.File(path, "r") as handle:
        return _h5_group_bytes(handle["/uns/spatial"]) if "/uns/spatial" in handle else 0


def _h5_image_count(path: str | Path) -> int:
    with h5py.File(path, "r") as handle:
        spatial = handle.get("/uns/spatial")
        if not isinstance(spatial, h5py.Group):
            return 0
        count = 0
        for library in spatial.values():
            if isinstance(library, h5py.Group) and isinstance(library.get("images"), h5py.Group):
                count += len(library["images"])
        return count


def _iter_matrix_blocks(matrix, n_obs: int, chunk_size: int = 2048) -> Iterable[Any]:
    if sparse.issparse(matrix):
        yield matrix
        return
    if isinstance(matrix, np.ndarray):
        yield matrix
        return
    for start in range(0, n_obs, chunk_size):
        yield matrix[start : min(start + chunk_size, n_obs)]


def _matrix_is_nonnegative_integer(matrix, n_obs: int) -> bool:
    for block in _iter_matrix_blocks(matrix, n_obs):
        values = block.data if sparse.issparse(block) else np.asarray(block).ravel()
        if values.size == 0:
            continue
        if not np.isfinite(values).all() or float(values.min()) < -1e-8:
            return False
        if not np.allclose(values, np.rint(values), rtol=0, atol=1e-6):
            return False
    return True


def _inspect_matrix(path: str | Path, requested: str) -> dict[str, Any]:
    backed = ad.read_h5ad(path, backed="r")
    try:
        source = requested
        if requested == "auto":
            source = "raw" if backed.raw is not None else "X"
        if source == "raw":
            if backed.raw is None:
                raise ValueError(f"seurat_matrix=raw was requested but {path} has no .raw")
            matrix = backed.raw.X
            n_vars = backed.raw.n_vars
            is_integer = _matrix_is_nonnegative_integer(matrix, backed.n_obs)
            if not is_integer:
                raise ValueError(f"{path}: .raw.X is not a valid non-negative integer count matrix")
            target_layer = "counts"
        else:
            matrix = backed.X
            n_vars = backed.n_vars
            is_integer = _matrix_is_nonnegative_integer(matrix, backed.n_obs)
            target_layer = "counts" if is_integer else "data"
        return {
            "input": str(path),
            "source": source,
            "target_layer": target_layer,
            "observations": backed.n_obs,
            "genes": n_vars,
            "matrix_bytes": _h5_matrix_bytes(path, source),
            "image_bytes": _h5_spatial_bytes(path),
            "images": _h5_image_count(path),
        }
    finally:
        backed.file.close()


def _resolve_seurat_selections(paths: Sequence[str], requested: str) -> list[dict[str, Any]]:
    normalized = requested.strip().lower()
    if normalized not in {"auto", "raw", "x"}:
        raise ValueError("seurat_matrix must be auto, raw, or X")
    source = "X" if normalized == "x" else normalized
    selections = [_inspect_matrix(path, source) for path in paths]
    layers = {selection["target_layer"] for selection in selections}
    if len(layers) != 1:
        details = ", ".join(
            f"{Path(item['input']).name}:{item['source']}->{item['target_layer']}"
            for item in selections
        )
        raise ValueError(
            "Cannot combine raw-count and normalized-only inputs into one Seurat assay. "
            f"Convert them separately or choose a consistent matrix source: {details}"
        )
    return selections


def _available_memory_bytes() -> int:
    candidates: list[int] = []
    for limit_path, current_path in (
        (Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current")),
        (Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"), Path("/sys/fs/cgroup/memory/memory.usage_in_bytes")),
    ):
        try:
            limit_text = limit_path.read_text().strip()
            if limit_text != "max":
                limit = int(limit_text)
                current = int(current_path.read_text().strip())
                if 0 < current < limit < 1 << 60:
                    candidates.append(limit - current)
        except (OSError, ValueError):
            pass
    try:
        pages = os.sysconf("SC_AVPHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        candidates.append(int(pages) * int(page_size))
    except (OSError, ValueError):
        pass
    return min(candidates) if candidates else 0


def _check_memory(selections: Sequence[dict[str, Any]], memory_limit_gb: float) -> tuple[int, int]:
    matrix_bytes = sum(int(item["matrix_bytes"]) for item in selections)
    image_bytes = sum(int(item["image_bytes"]) for item in selections)
    estimated = 5 * matrix_bytes + image_bytes + GIB
    available = _available_memory_bytes()
    if memory_limit_gb > 0:
        allowed = int(memory_limit_gb * GIB)
    elif available > 0:
        allowed = int(available * 0.8)
    else:
        allowed = 0
    logger.info(
        "Seurat memory preflight: matrix=%.2f GiB images=%.2f GiB estimated_peak=%.2f GiB allowed=%s",
        matrix_bytes / GIB,
        image_bytes / GIB,
        estimated / GIB,
        f"{allowed / GIB:.2f} GiB" if allowed else "unknown",
    )
    if allowed and estimated > allowed:
        raise MemoryError(
            f"Estimated Seurat conversion peak ({estimated / GIB:.2f} GiB) exceeds the "
            f"allowed memory ({allowed / GIB:.2f} GiB). Disable images, split samples, "
            "or increase --memory_limit_gb."
        )
    return estimated, allowed


def _materialize_matrix(matrix, n_obs: int):
    if sparse.issparse(matrix):
        return matrix.tocsr(copy=True)
    if isinstance(matrix, np.ndarray):
        return np.array(matrix, copy=True)
    blocks = []
    for block in _iter_matrix_blocks(matrix, n_obs):
        blocks.append(block.tocsr() if sparse.issparse(block) else np.asarray(block))
    if any(sparse.issparse(block) for block in blocks):
        return sparse.vstack(
            [block if sparse.issparse(block) else sparse.csr_matrix(block) for block in blocks],
            format="csr",
        )
    return np.concatenate(blocks, axis=0)


def _minimal_seurat_h5ad(
    source_path: str,
    selection: dict[str, Any],
    output: Path,
    used_library_ids: set[str],
) -> list[str]:
    backed = ad.read_h5ad(source_path, backed="r")
    try:
        source = backed.raw if selection["source"] == "raw" else backed
        matrix = _materialize_matrix(source.X, backed.n_obs)
        minimal = ad.AnnData(
            X=matrix,
            obs=backed.obs.copy(),
            var=source.var.copy(),
        )
        for key, value in backed.obsm.items():
            try:
                minimal.obsm[key] = value.to_memory() if hasattr(value, "to_memory") else np.asarray(value)
            except (TypeError, ValueError):
                logger.warning("Skipping unsupported obsm[%s] during Seurat export", key)
        if "spatial" in backed.uns:
            minimal.uns["spatial"] = copy.deepcopy(backed.uns["spatial"])
        source_name = _clean_stem(source_path)
        original_libraries = (
            list(dict.fromkeys(minimal.obs["library_id"].astype(str)))
            if "library_id" in minimal.obs
            else [source_name]
        )
        mapping = {
            library: _unique_library_id(library, source_name, used_library_ids)
            for library in original_libraries
        }
        if "library_id" in minimal.obs:
            minimal.obs["library_id"] = pd.Categorical(
                minimal.obs["library_id"].astype(str).map(mapping)
            )
        else:
            minimal.obs["library_id"] = pd.Categorical([mapping[source_name]] * minimal.n_obs)

        spatial = minimal.uns.get("spatial")
        if isinstance(spatial, dict) and spatial:
            if len(spatial) == 1 and len(mapping) == 1:
                minimal.uns["spatial"] = {next(iter(mapping.values())): next(iter(spatial.values()))}
            else:
                renamed_spatial = {}
                for key, entry in spatial.items():
                    target = mapping.get(str(key))
                    if target is None:
                        target = _unique_library_id(str(key), source_name, used_library_ids)
                    renamed_spatial[target] = entry
                minimal.uns["spatial"] = renamed_spatial
        minimal.write_h5ad(output)
        return list(mapping.values())
    finally:
        backed.file.close()


def _native_spatial_h5ad(path: str | Path) -> bool:
    backed = ad.read_h5ad(path, backed="r")
    try:
        if "spatial" not in backed.obsm or "library_id" not in backed.obs:
            return False
        spatial = backed.uns.get("spatial")
        if not isinstance(spatial, dict) or not spatial:
            return False
        libraries = set(backed.obs["library_id"].astype(str))
        for library in libraries:
            entry = spatial.get(library)
            if not isinstance(entry, dict):
                return False
            if not entry.get("images") or not entry.get("scalefactors"):
                return False
        return True
    finally:
        backed.file.close()


def _h5ad_has_spatial_coordinates(path: str | Path) -> bool:
    backed = ad.read_h5ad(path, backed="r")
    try:
        return "spatial" in backed.obsm
    finally:
        backed.file.close()


def _run_schard(
    input_h5ad: Path,
    output_rds: Path,
    data_type: str,
    target_layer: str,
    conversion_mode: str,
) -> None:
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "h5ad_to_seurat.R"
    staged = output_rds.parent / f".{output_rds.name}.tmp-{os.getpid()}"
    if staged.exists():
        staged.unlink()
    command = [
        "Rscript",
        str(script_path),
        str(input_h5ad),
        str(staged),
        data_type,
        target_layer,
        conversion_mode,
        SCHARD_MIN_VERSION,
    ]
    logger.info("Running Schard: %s", " ".join(command))
    subprocess.run(command, check=True)
    if not staged.exists() or staged.stat().st_size == 0:
        raise RuntimeError("Schard completed without creating a valid RDS file")
    _replace_file(staged, output_rds)


def h5ad_to_seurat(
    INPUT_list,
    output_dir,
    data_type="auto",
    seurat_matrix="auto",
    memory_limit_gb=0,
    output_name="",
):
    paths = [str(value) for value in INPUT_list]
    if not paths:
        raise ValueError("At least one H5AD input is required")
    selections = _resolve_seurat_selections(paths, seurat_matrix)
    estimated, allowed = _check_memory(selections, float(memory_limit_gb))
    target_layer = selections[0]["target_layer"]
    output_directory = Path(output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    base_name = output_name or _clean_stem(paths[0])
    output_rds = output_directory / f"{base_name}.rds"

    with tempfile.TemporaryDirectory(prefix="spatialsnake-seurat-", dir=output_directory) as temp_dir:
        temp_paths: list[Path] = []
        used_library_ids: set[str] = set()
        path_suffixes: list[str] = []
        used_suffixes: set[str] = set()
        for index, (path, selection) in enumerate(zip(paths, selections, strict=True)):
            temp_path = Path(temp_dir) / f"selected-{index}.h5ad"
            _minimal_seurat_h5ad(path, selection, temp_path, used_library_ids)
            temp_paths.append(temp_path)
            path_suffixes.append(_unique_library_id(_clean_stem(path), "input", used_suffixes))
        _resolve_duplicate_obs_names(temp_paths, path_suffixes)
        combined = Path(temp_dir) / "seurat_input.h5ad"
        _concat_h5ad_on_disk(temp_paths, combined)
        combined_summary = _validate_h5ad(
            combined, sum(item["observations"] for item in selections)
        )

        normalized_type = data_type.strip().lower()
        if normalized_type == "auto":
            normalized_type = "st" if _h5ad_has_spatial_coordinates(combined) else "sc"
        if normalized_type not in {"st", "sc"}:
            raise ValueError("type must be auto, st, or sc")
        native = normalized_type == "st" and _native_spatial_h5ad(combined)
        conversion_mode = "spatial" if native else "generic"
        conversion_warning = ""
        if normalized_type == "st" and not native:
            conversion_warning = (
                "complete Visium-style image metadata was unavailable; used generic Seurat "
                "conversion and retained spatial coordinates as a dimensional reduction"
            )
            logger.warning(
                "No complete Visium-style image metadata was found; exporting a generic Seurat "
                "object with spatial coordinates as a dimensional reduction."
            )
        _run_schard(
            combined,
            output_rds,
            normalized_type,
            target_layer,
            conversion_mode,
        )

    report = pd.DataFrame(selections)
    report["output"] = str(output_rds)
    report["data_type"] = normalized_type
    report["conversion_mode"] = conversion_mode
    report["libraries"] = len(combined_summary["libraries"]) or len(paths)
    report["estimated_peak_gib"] = estimated / GIB
    report["allowed_memory_gib"] = allowed / GIB if allowed else np.nan
    report["warnings"] = conversion_warning
    return output_rds, report


def h5ad_to_zarr(INPUT_list, output_zarr, base_name=None):
    paths = [str(value) for value in INPUT_list]
    if len(paths) != 1:
        raise ValueError(
            "h5ad -> zarr accepts exactly one input. Use useful_tool --option=merge "
            "before or after conversion when multiple objects are required."
        )
    input_path = paths[0]
    adata = ad.read_h5ad(input_path)
    if "spatial" in adata.obsm and isinstance(adata.uns.get("spatial"), dict):
        sdata = from_legacy_anndata(adata)
        mode = "legacy_spatial"
    else:
        sdata = sd.SpatialData(tables={"table": adata.copy()})
        mode = "table_only"

    output = Path(output_zarr)
    output.parent.mkdir(parents=True, exist_ok=True)
    staged = output.parent / f".{output.name}.tmp-{os.getpid()}"
    if staged.exists():
        shutil.rmtree(staged)
    sdata.write(staged)
    reopened = sd.read_zarr(staged)
    _, table = _select_table(reopened, "")
    if table.n_obs != adata.n_obs or not table.obs_names.astype(str).equals(adata.obs_names.astype(str)):
        raise ValueError("H5AD -> Zarr validation failed: observation IDs or order changed")
    if output.exists():
        shutil.rmtree(output)
    os.replace(staged, output)
    return sdata, pd.DataFrame(
        [{
            "input": input_path,
            "output": str(output),
            "observations": adata.n_obs,
            "genes": adata.n_vars,
            "conversion_mode": mode,
            "matrix_source": "X+raw+layers",
            "libraries": 1,
            "warnings": "",
        }]
    )


def _validate_formats(inputs: Sequence[str], transform_from: str, transform_to: str) -> None:
    source = transform_from.lower()
    target = transform_to.lower()
    allowed = {("zarr", "h5ad"), ("zarr", "seurat"), ("h5ad", "zarr"), ("h5ad", "seurat")}
    if (source, target) not in allowed:
        raise ValueError(
            f"Unsupported conversion {source!r} -> {target!r}. Allowed conversions: "
            "zarr->h5ad, zarr->seurat, h5ad->zarr, h5ad->seurat."
        )
    for path in inputs:
        lower = str(path).lower().rstrip("/")
        if source == "h5ad" and not lower.endswith(".h5ad"):
            raise ValueError(f"Expected an .h5ad input, got {path}")
        if source == "zarr" and not Path(path).is_dir():
            raise ValueError(f"Expected a SpatialData Zarr directory, got {path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Spatialsnake format conversion")
    parser.add_argument("--INPUT", nargs="+", required=True, help="Input Zarr or H5AD path(s)")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--save_image", default="True")
    parser.add_argument("--transform_to", required=True)
    parser.add_argument("--transform_from", required=True)
    parser.add_argument("--type", default="auto", choices=["auto", "st", "sc"])
    parser.add_argument("--seurat_matrix", default="auto")
    parser.add_argument("--table_key", default="")
    parser.add_argument("--memory_limit_gb", type=float, default=0)
    parser.add_argument("--keep_intermediate", default="False")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    inputs = [str(value) for value in args.INPUT]
    source = args.transform_from.strip().lower()
    target = args.transform_to.strip().lower()
    _validate_formats(inputs, source, target)
    save_image = parse_bool(args.save_image)
    keep_intermediate = parse_bool(args.keep_intermediate)
    if args.memory_limit_gb < 0:
        raise ValueError("memory_limit_gb must be >= 0")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base_name = _clean_stem(inputs[0])
    report: pd.DataFrame
    log_step(logger, 1, 3, f"loading {len(inputs)} {source} input file(s)")

    if source == "zarr" and target == "h5ad":
        output = output_dir / f"{base_name}.h5ad"
        log_step(logger, 2, 3, f"exporting SpatialData to H5AD: {output}")
        report = _export_zarr_to_h5ad(inputs, output, save_image, args.table_key)
        report["output"] = str(output)
        report["matrix_source"] = "X+raw+layers"
        report["conversion_mode"] = "legacy_h5ad"
        report["libraries"] = len(report)
    elif source == "zarr" and target == "seurat":
        intermediate = output_dir / f"{base_name}.h5ad" if keep_intermediate else None
        with tempfile.TemporaryDirectory(prefix="spatialsnake-transform-rds-", dir=output_dir) as temp_dir:
            h5ad_path = intermediate or Path(temp_dir) / f"{base_name}.h5ad"
            log_step(logger, 2, 3, "exporting minimal-compatible H5AD and converting with Schard")
            zarr_report = _export_zarr_to_h5ad(inputs, h5ad_path, save_image, args.table_key)
            output, report = h5ad_to_seurat(
                [str(h5ad_path)],
                output_dir,
                data_type=args.type,
                seurat_matrix=args.seurat_matrix,
                memory_limit_gb=args.memory_limit_gb,
                output_name=base_name,
            )
            report["zarr_libraries"] = len(zarr_report)
            report["input"] = ";".join(inputs)
            report["intermediate_h5ad"] = str(intermediate) if intermediate else ""
    elif source == "h5ad" and target == "zarr":
        output = output_dir / f"{base_name}.zarr"
        log_step(logger, 2, 3, f"converting H5AD to SpatialData Zarr: {output}")
        _, report = h5ad_to_zarr(inputs, output, base_name)
    else:
        log_step(logger, 2, 3, "converting H5AD to Seurat with Schard")
        output, report = h5ad_to_seurat(
            inputs,
            output_dir,
            data_type=args.type,
            seurat_matrix=args.seurat_matrix,
            memory_limit_gb=args.memory_limit_gb,
            output_name=base_name,
        )

    report_path = output_dir / f"{base_name}_transform_report.csv"
    report.to_csv(report_path, index=False)
    log_step(logger, 3, 3, f"validated output and wrote report: {report_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        logger.error("Transform failed: %s", exc)
        raise
