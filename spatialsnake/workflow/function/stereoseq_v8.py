from __future__ import annotations

import os
import re
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any

import anndata as ad
import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
from shapely.geometry import Point, Polygon
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Identity

try:
    from dask_image.imread import imread
except ImportError:
    try:
        import tifffile

        def imread(path, **kwargs):
            return tifffile.imread(path)

    except ImportError:
        from PIL import Image

        def imread(path, **kwargs):
            return np.asarray(Image.open(path))

__all__ = ["stereoseq_v8"]


class SK:
    GENE_EXP = "geneExp"
    EXPRESSION = "expression"
    FEATURE_KEY = "gene"
    GENE_NAME = "geneName"
    OFFSET = "offset"
    COUNT = "count"
    COORD_X = "x"
    COORD_Y = "y"
    RESOLUTION = "resolution"
    REGION_KEY = "region"
    INSTANCE_KEY = "instance_id"


def _decode_value(x: Any) -> Any:
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode("utf-8", errors="ignore")
    return x


def _decode_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if out[col].dtype == object or getattr(out[col].dtype, "kind", None) == "S":
            out[col] = out[col].map(_decode_value)
    return out


def _to_dataframe(obj: h5py.Dataset | np.ndarray, fallback_columns: list[str] | None = None) -> pd.DataFrame:
    arr = obj[:] if isinstance(obj, h5py.Dataset) else obj
    if getattr(arr.dtype, "fields", None) is not None:
        df = pd.DataFrame.from_records(arr)
    else:
        df = pd.DataFrame(arr, columns=fallback_columns)
    return _decode_dataframe(df)


def _ensure_columns(df: pd.DataFrame, expected: list[str]) -> pd.DataFrame:
    out = df.copy()
    missing = [c for c in expected if c not in out.columns]
    if not missing:
        return out
    if len(out.columns) < len(expected):
        raise ValueError(f"Expected columns {expected}, got {list(out.columns)}")
    rename_map = {}
    for old, new in zip(list(out.columns)[: len(expected)], expected):
        if new not in out.columns:
            rename_map[old] = new
    return out.rename(columns=rename_map)


def _make_unique_gene_names(gene_names: list[str], gene_ids: list[str]) -> list[str]:
    counts = pd.Series(gene_names).value_counts()
    out = []
    for gname, gid in zip(gene_names, gene_ids):
        out.append(f"{gname}_{gid}" if counts[gname] > 1 else gname)
    return out


def _normalize_bin_requests(bin_sizes: list[int | str] | int | str | None) -> tuple[list[int] | None, str | None]:
    if bin_sizes is None:
        return None, None
    items = [bin_sizes] if isinstance(bin_sizes, (int, str)) else list(bin_sizes)
    tissue_bins: list[int] = []
    cellbin_mode: str | None = None
    for item in items:
        if isinstance(item, str):
            key = item.strip().lower()
            if key in {"cellbin", "cell_bin", "cell"}:
                cellbin_mode = "cellbin"
                continue
            if key in {"adjusted_cellbin", "adjusted.cellbin", "adjusted-cellbin", "adjusted"}:
                cellbin_mode = "adjusted_cellbin"
                continue
        tissue_bins.append(int(item))
    return tissue_bins, cellbin_mode


def _infer_dataset_id(feature_exp_path: Path) -> str:
    for pat in ["*.adjusted.cellbin.gef", "*.cellbin.gef", "*.tissue.gef", "*.gef"]:
        hits = sorted(feature_exp_path.glob(pat))
        if hits:
            name = hits[0].name
            for suffix in [".adjusted.cellbin.gef", ".cellbin.gef", ".tissue.gef", ".gef"]:
                if name.endswith(suffix):
                    return name[: -len(suffix)]
    raise ValueError(f"Could not infer dataset_id from files under {feature_exp_path}")


def _find_cellbin_path(feature_exp_path: Path, dataset_id: str, adjusted: bool = False) -> Path:
    candidates = (
        [
            feature_exp_path / f"{dataset_id}.adjusted.cellbin.gef",
            feature_exp_path / f"{dataset_id}.cellbin.gef",
        ]
        if adjusted
        else [
            feature_exp_path / f"{dataset_id}.cellbin.gef",
            feature_exp_path / f"{dataset_id}.adjusted.cellbin.gef",
        ]
    )
    for path in candidates:
        if path.exists():
            return path
    glob_candidates = sorted(feature_exp_path.glob("*cellbin*.gef"))
    if not glob_candidates:
        raise ValueError(f"No cellbin.gef file found under {feature_exp_path}")
    return glob_candidates[0]


def _default_transformations() -> dict[str, Identity]:
    return {"global": Identity()}


def _read_images(
    image_dir: Path,
    imread_kwargs: Mapping[str, Any],
    image_models_kwargs: Mapping[str, Any],
) -> dict[str, Any]:
    if not image_dir.exists():
        return {}
    parse_kwargs = dict(image_models_kwargs)
    parse_kwargs.setdefault("transformations", _default_transformations())
    patterns = [re.compile(r".*_HE_regist\.tif$"), re.compile(r".*_HE_tissue_cut\.tif$")]
    images = {}
    for image_file in os.listdir(image_dir):
        if not any(pattern.match(image_file) for pattern in patterns):
            continue
        image_name = Path(image_file).stem
        img_data = imread(image_dir / image_file, **imread_kwargs)
        if img_data.ndim == 4 and img_data.shape[0] == 1:
            img_data = img_data[0]
        if img_data.ndim == 3 and img_data.shape[-1] in [1, 3, 4]:
            img_data = np.moveaxis(img_data, -1, 0)
        images[image_name] = Image2DModel.parse(img_data, dims=("c", "y", "x"), **parse_kwargs)
    return images


def _read_tissue_bins(
    tissue_gef_path: Path,
    requested_bin_sizes: list[int] | None,
    tables: dict[str, Any],
    points: dict[str, Any],
) -> None:
    with h5py.File(str(tissue_gef_path), "r") as tissue_gef:
        available_bins = list(tissue_gef[SK.GENE_EXP].keys())
        bins_to_read = available_bins if requested_bin_sizes is None else [f"bin{size}" for size in requested_bin_sizes]
        missing_bins = set(bins_to_read) - set(available_bins)
        if missing_bins:
            raise ValueError(f"Requested bin sizes {missing_bins} not found in tissue.gef. Available bins: {available_bins}")

        for bin_name in bins_to_read:
            bin_group = tissue_gef[SK.GENE_EXP][bin_name]
            expr_ds = bin_group[SK.EXPRESSION]
            resolution = dict(expr_ds.attrs).get(SK.RESOLUTION)
            if isinstance(resolution, np.ndarray):
                resolution = resolution.item()

            gene_df = _ensure_columns(_to_dataframe(bin_group[SK.FEATURE_KEY], ["geneID", SK.GENE_NAME, SK.OFFSET, SK.COUNT]), ["geneID", SK.GENE_NAME, SK.OFFSET, SK.COUNT])
            expr_df = _ensure_columns(_to_dataframe(expr_ds, [SK.COORD_X, SK.COORD_Y, SK.COUNT]), [SK.COORD_X, SK.COORD_Y, SK.COUNT])

            gene_df["geneID"] = gene_df["geneID"].astype(str)
            gene_df[SK.GENE_NAME] = gene_df[SK.GENE_NAME].astype(str)
            gene_df["gene_name_unique"] = _make_unique_gene_names(gene_df[SK.GENE_NAME].tolist(), gene_df["geneID"].tolist())

            offset_ok = np.array_equal(
                gene_df[SK.OFFSET].to_numpy(),
                np.insert(np.cumsum(gene_df[SK.COUNT].to_numpy()), 0, 0)[:-1],
            )
            if not offset_ok:
                raise ValueError(f"Offset validation failed for {tissue_gef_path} / {bin_name}")

            repeated_gene_names = np.repeat(gene_df["gene_name_unique"].to_numpy(), gene_df[SK.COUNT].to_numpy(dtype=np.int64))
            if len(repeated_gene_names) != len(expr_df):
                raise ValueError(f"Expression length mismatch in {bin_name}: {len(repeated_gene_names)} vs {len(expr_df)}")

            expr_df[SK.FEATURE_KEY] = repeated_gene_names
            gene_codes, genes_with_expression = pd.factorize(expr_df[SK.FEATURE_KEY], sort=False)
            points_coords = expr_df[[SK.COORD_X, SK.COORD_Y]].drop_duplicates().reset_index(drop=True)
            points_coords["bin_id"] = np.arange(len(points_coords), dtype=np.int64)
            index_to_bin_id = pd.merge(
                expr_df[[SK.COORD_X, SK.COORD_Y]],
                points_coords,
                on=[SK.COORD_X, SK.COORD_Y],
                how="left",
                validate="many_to_one",
            )
            expression_matrix = coo_matrix(
                (
                    expr_df[SK.COUNT].to_numpy(dtype=np.float32),
                    (index_to_bin_id["bin_id"].to_numpy(dtype=np.int32), gene_codes.astype(np.int32)),
                ),
                shape=(len(points_coords), len(genes_with_expression)),
            ).tocsr()

            point_ids = pd.Index([str(i) for i in range(len(points_coords))], name=None)
            points_df = points_coords.drop(columns=["bin_id"]).copy()
            points_df.index = point_ids
            points_key = f"{bin_name}_genes"
            table_key = f"{bin_name}_table"

            obs = pd.DataFrame({SK.INSTANCE_KEY: point_ids, SK.REGION_KEY: pd.Categorical([points_key] * len(point_ids))}, index=point_ids)
            var = gene_df.set_index("gene_name_unique").loc[list(genes_with_expression), :].copy()
            var.index.name = None
            adata = ad.AnnData(expression_matrix, obs=obs, var=var)
            adata.uns["resolution"] = resolution
            adata.uns["bin_name"] = bin_name

            tables[table_key] = TableModel.parse(adata, region=points_key, region_key=SK.REGION_KEY, instance_key=SK.INSTANCE_KEY)
            points[points_key] = PointsModel.parse(points_df, coordinates={"x": SK.COORD_X, "y": SK.COORD_Y}, transformations=_default_transformations())


def _clean_cellbin_border_offsets(offsets: np.ndarray, sentinel: float = 32767.0, max_abs_offset: float = 256.0) -> np.ndarray:
    offsets = np.asarray(offsets, dtype=np.float32)
    valid = np.isfinite(offsets).all(axis=1)
    valid &= ~(np.isclose(offsets[:, 0], sentinel) | np.isclose(offsets[:, 1], sentinel))
    valid &= np.abs(offsets[:, 0]) <= max_abs_offset
    valid &= np.abs(offsets[:, 1]) <= max_abs_offset
    cleaned = offsets[valid]
    if len(cleaned) == 0:
        return cleaned
    dedup = [tuple(cleaned[0])]
    for xy in cleaned[1:]:
        point = tuple(xy)
        if point != dedup[-1]:
            dedup.append(point)
    cleaned = np.asarray(dedup, dtype=np.float32)
    if len(cleaned) >= 2 and np.allclose(cleaned[0], cleaned[-1]):
        cleaned = cleaned[:-1]
    return cleaned


def _build_cellbin_shapes(cell_df: pd.DataFrame, cell_ids: pd.Index, cell_border: np.ndarray | None) -> gpd.GeoDataFrame:
    geometries = []
    if cell_border is not None and len(cell_border) != len(cell_df):
        raise ValueError(f"cellBorder length mismatch: len(cellBorder)={len(cell_border)} vs len(cell_df)={len(cell_df)}")
    for i, (_, row) in enumerate(cell_df.iterrows()):
        cx = float(row["x"])
        cy = float(row["y"])
        area = float(row["area"]) if "area" in row and pd.notna(row["area"]) else 16.0
        poly = None
        if cell_border is not None:
            offsets = _clean_cellbin_border_offsets(cell_border[i])
            if len(offsets) >= 3:
                coords = np.column_stack([cx + offsets[:, 0], cy + offsets[:, 1]])
                valid_abs = np.isfinite(coords).all(axis=1)
                valid_abs &= coords[:, 0] > 0
                valid_abs &= coords[:, 1] > 0
                valid_abs &= coords[:, 0] < 50000
                valid_abs &= coords[:, 1] < 50000
                coords = coords[valid_abs]
                if len(coords) >= 3:
                    dedup = []
                    for xy in coords.tolist():
                        point = tuple(xy)
                        if not dedup or point != dedup[-1]:
                            dedup.append(point)
                    if len(set(dedup)) >= 3:
                        try:
                            poly = Polygon(dedup)
                            if not poly.is_valid:
                                poly = poly.buffer(0)
                            if poly.is_empty or poly.area <= 0:
                                poly = None
                        except Exception:
                            poly = None
        if poly is None:
            radius = max(np.sqrt(max(area, 1.0) / np.pi), 1.0)
            poly = Point(cx, cy).buffer(radius, resolution=16)
        geometries.append(poly)
    return gpd.GeoDataFrame({"cell_id": cell_ids, "geometry": geometries}, index=cell_ids)


def _map_cell_type_names(cell_df: pd.DataFrame, cell_type_list: np.ndarray | None) -> pd.DataFrame:
    out = cell_df.copy()
    if cell_type_list is None or "cellTypeID" not in out.columns:
        return out
    labels = [_decode_value(x) for x in np.asarray(cell_type_list)]
    mapping = {i: str(name) for i, name in enumerate(labels)}
    out["cellTypeName"] = out["cellTypeID"].map(mapping).fillna(out["cellTypeID"].astype(str))
    out["cellTypeName"] = out["cellTypeName"].astype("category")
    return out


def _read_cellbin(cellbin_gef_path: Path, tables: dict[str, Any], shapes: dict[str, Any], points: dict[str, Any]) -> None:
    mode_prefix = "adjusted_cellbin" if "adjusted.cellbin" in cellbin_gef_path.name else "cellbin"
    table_key = f"{mode_prefix}_table"
    shape_key = f"{mode_prefix}_shapes"
    point_key = f"{mode_prefix}_points"
    with h5py.File(str(cellbin_gef_path), "r") as handle:
        if "cellBin" not in handle:
            raise ValueError(f"'cellBin' group not found in {cellbin_gef_path}")
        group = handle["cellBin"]
        if "cell" not in group or "gene" not in group or "cellExp" not in group:
            raise ValueError(f"{cellbin_gef_path} is missing one of required datasets: 'cell', 'gene', 'cellExp'")
        cell_df = _to_dataframe(group["cell"])
        gene_df = _to_dataframe(group["gene"])
        cell_exp_df = _to_dataframe(group["cellExp"])
        cell_border = np.asarray(group["cellBorder"][:]) if "cellBorder" in group else None
        cell_type_list = np.asarray(group["cellTypeList"][:]) if "cellTypeList" in group else None

    cell_df = _map_cell_type_names(_ensure_columns(cell_df, ["id", "x", "y", "offset", "geneCount", "expCount", "dnbCount", "area", "cellTypeID", "clusterID"]), cell_type_list)
    gene_df = _ensure_columns(gene_df, ["geneID", "geneName", "offset", "cellCount", "expCount", "maxMIDcount"])
    cell_exp_df = _ensure_columns(cell_exp_df, ["geneID", "count"])

    gene_df["geneID"] = gene_df["geneID"].astype(str)
    gene_df["geneName"] = gene_df["geneName"].astype(str)
    gene_df["gene_name_unique"] = _make_unique_gene_names(gene_df["geneName"].tolist(), gene_df["geneID"].tolist())
    gene_df.index = pd.Index(gene_df["gene_name_unique"].tolist(), name=None)

    offsets = cell_df["offset"].to_numpy(dtype=np.int64, copy=False)
    gene_count = cell_df["geneCount"].to_numpy(dtype=np.int64, copy=False)
    total_nnz = int(np.sum(gene_count))
    raw_gene_ids = cell_exp_df["geneID"].to_numpy(dtype=np.int64, copy=False)
    raw_counts = cell_exp_df["count"].to_numpy(dtype=np.float32, copy=False)
    if total_nnz <= 0:
        raise ValueError(f"geneCount sums to zero in {cellbin_gef_path}")
    if total_nnz > len(raw_gene_ids):
        raise ValueError(f"cellExp is shorter than expected: total_nnz={total_nnz}, len(cellExp)={len(raw_gene_ids)}")

    n_genes = len(gene_df)
    gid_max = int(np.max(raw_gene_ids[:total_nnz]))
    gid_min = int(np.min(raw_gene_ids[:total_nnz]))
    adjust = 1 if gid_max == n_genes or (gid_max != n_genes - 1 and gid_min == 1 and 0 not in raw_gene_ids[: min(1000, len(raw_gene_ids))]) else 0
    indices = raw_gene_ids[:total_nnz] - adjust
    if int(np.min(indices)) < 0 or int(np.max(indices)) >= n_genes:
        raise ValueError(f"Mapped geneID out of range in {cellbin_gef_path}: min={int(np.min(indices))}, max={int(np.max(indices))}, n_genes={n_genes}")

    expected_offsets = np.insert(np.cumsum(gene_count), 0, 0)[:-1]
    if not np.array_equal(offsets, expected_offsets):
        raise ValueError(f"Offset validation failed for {cellbin_gef_path}. This reader assumes row-wise contiguous cellExp slices.")

    indptr = np.insert(np.cumsum(gene_count), 0, 0).astype(np.int64)
    X = csr_matrix((raw_counts[:total_nnz], indices.astype(np.int32), indptr), shape=(len(cell_df), n_genes))
    cell_ids = pd.Index(cell_df["id"].astype(str).tolist(), name=None)
    obs = cell_df.copy()
    obs.index = cell_ids
    obs["cell_id"] = cell_ids
    obs[SK.INSTANCE_KEY] = cell_ids
    obs[SK.REGION_KEY] = pd.Categorical([shape_key] * len(cell_ids))
    for col in ["clusterID", "cellTypeID", "cellTypeName"]:
        if col in obs.columns:
            obs[col] = obs[col].astype("category")

    adata = ad.AnnData(X=X, obs=obs.copy(), var=gene_df.copy())
    adata.obsm["spatial"] = obs[["x", "y"]].to_numpy(dtype=np.float32)
    adata.uns["source_gef"] = str(cellbin_gef_path)
    adata.uns["mode"] = mode_prefix

    shapes[shape_key] = ShapesModel.parse(_build_cellbin_shapes(obs, cell_ids, cell_border), transformations=_default_transformations())
    tables[table_key] = TableModel.parse(adata, region=shape_key, region_key=SK.REGION_KEY, instance_key=SK.INSTANCE_KEY)

    point_df = obs[["x", "y"]].copy()
    point_df.index = cell_ids
    points[point_key] = PointsModel.parse(point_df, coordinates={"x": "x", "y": "y"}, transformations=_default_transformations())


def _load_analysis_tables(analysis_dir: Path, tables: dict[str, Any], points: dict[str, Any]) -> None:
    if not analysis_dir.exists():
        return
    for h5ad_file in sorted(analysis_dir.glob("*.h5ad")):
        adata_analysis = ad.read_h5ad(h5ad_file)
        if "spatial" not in adata_analysis.obsm:
            continue
        adata_analysis.uns.pop("spatialdata_attrs", None)
        region_name = f"analysis_{h5ad_file.stem}_points"
        adata_analysis.obs_names = adata_analysis.obs_names.astype(str)
        adata_analysis.obs[SK.INSTANCE_KEY] = adata_analysis.obs_names
        adata_analysis.obs[SK.REGION_KEY] = pd.Categorical([region_name] * adata_analysis.n_obs)
        tables[f"analysis_{h5ad_file.stem}"] = TableModel.parse(
            adata_analysis,
            region=region_name,
            region_key=SK.REGION_KEY,
            instance_key=SK.INSTANCE_KEY,
        )
        points_df = pd.DataFrame(adata_analysis.obsm["spatial"], columns=[SK.COORD_X, SK.COORD_Y], index=adata_analysis.obs_names)
        points[region_name] = PointsModel.parse(points_df, coordinates={"x": SK.COORD_X, "y": SK.COORD_Y}, transformations=_default_transformations())


def stereoseq_v8(
    path: str | Path,
    dataset_id: str | None = None,
    bin_sizes: list[int | str] | int | str | None = None,
    load_analysis: bool = True,
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> SpatialData:
    path = Path(path)
    feature_exp_path = path / "feature_expression"
    image_dir = path / "image"
    analysis_dir = path / "analysis"
    if not feature_exp_path.exists():
        raise ValueError(f"feature_expression directory not found in {path}")
    if dataset_id is None:
        dataset_id = _infer_dataset_id(feature_exp_path)

    requested_tissue_bins, cellbin_mode = _normalize_bin_requests(bin_sizes)
    images = _read_images(image_dir=image_dir, imread_kwargs=imread_kwargs, image_models_kwargs=image_models_kwargs)
    tables: dict[str, Any] = {}
    points: dict[str, Any] = {}
    shapes: dict[str, Any] = {}

    tissue_gef_path = feature_exp_path / f"{dataset_id}.tissue.gef"
    if requested_tissue_bins is None and cellbin_mode is None:
        if tissue_gef_path.exists():
            _read_tissue_bins(tissue_gef_path, requested_tissue_bins, tables, points)
        else:
            raise ValueError(f"No default tissue.gef found at {tissue_gef_path}. Pass bin_sizes='cellbin' if you want to read cellbin.gef only.")
    else:
        if requested_tissue_bins:
            if not tissue_gef_path.exists():
                raise ValueError(f"tissue.gef file not found: {tissue_gef_path}")
            _read_tissue_bins(tissue_gef_path, requested_tissue_bins, tables, points)
        if cellbin_mode is not None:
            _read_cellbin(
                cellbin_gef_path=_find_cellbin_path(feature_exp_path, dataset_id, adjusted=(cellbin_mode == "adjusted_cellbin")),
                tables=tables,
                shapes=shapes,
                points=points,
            )

    if load_analysis:
        _load_analysis_tables(analysis_dir=analysis_dir, tables=tables, points=points)

    return SpatialData(images=images, tables=tables, points=points, shapes=shapes)
