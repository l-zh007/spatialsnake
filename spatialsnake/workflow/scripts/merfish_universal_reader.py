from __future__ import annotations

import glob
import logging
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd
from spatialdata.models import PointsModel

try:
    import pyarrow.parquet as pq
except Exception:
    pq = None


log = logging.getLogger("merfish_universal_reader")


@dataclass
class MerfishPointsConfig:
    points_key: str = "transcripts"
    max_points: Optional[int] = None
    auto_max_points: int = 3_000_000
    auto_trigger_rows: int = 20_000_000
    auto_trigger_bytes: int = 2_000_000_000
    keep_cell_id: bool = False
    keep_fov: bool = False
    keep_qv: bool = False


def _dedup(values: List[str]) -> List[str]:
    return list(dict.fromkeys(values))


def _pick_largest_file(paths: List[str]) -> Optional[str]:
    files = [p for p in paths if os.path.isfile(p)]
    if not files:
        return None
    files.sort(key=lambda p: os.path.getsize(p), reverse=True)
    return files[0]


def find_transcripts_file(sample_dir: str) -> Optional[str]:
    patterns = [
        "detected_transcripts.parquet",
        "detected_transcripts*.parquet",
        "*detected_transcripts*.parquet",
        "*transcripts*.parquet",
        "detected_transcripts.csv",
        "detected_transcripts*.csv",
        "*detected_transcripts*.csv",
        "*transcripts*.csv",
    ]
    candidates: List[str] = []
    for pattern in patterns:
        candidates.extend(glob.glob(os.path.join(sample_dir, pattern)))
    return _pick_largest_file(_dedup(candidates))


def _pick_col(columns: List[str], aliases: List[str]) -> Optional[str]:
    for alias in aliases:
        if alias in columns:
            return alias
    lower = {c.lower(): c for c in columns}
    for alias in aliases:
        key = alias.lower()
        if key in lower:
            return lower[key]
    return None


def infer_transcript_columns(columns: List[str]) -> Dict[str, Optional[str]]:
    return {
        "x": _pick_col(columns, ["global_x", "x", "x_global", "center_x", "x_um", "xc", "xcoord", "x_coordinate"]),
        "y": _pick_col(columns, ["global_y", "y", "y_global", "center_y", "y_um", "yc", "ycoord", "y_coordinate"]),
        "z": _pick_col(columns, ["global_z", "z", "z_index", "z_global", "zc", "zcoord", "z_coordinate"]),
        "gene": _pick_col(columns, ["gene", "gene_name", "gene_symbol", "target", "target_name", "feature_name", "features"]),
        "cell_id": _pick_col(columns, ["cell_id", "cell", "cellid", "cellID", "cell_index", "cell_idx"]),
        "qv": _pick_col(columns, ["qv", "quality", "quality_score", "confidence", "conf"]),
        "fov": _pick_col(columns, ["fov", "fov_id", "tile", "tile_id", "fovindex"]),
    }


def _auto_cap(path: str, cfg: MerfishPointsConfig) -> Optional[int]:
    if cfg.max_points is not None:
        return cfg.max_points
    size = 0
    try:
        size = os.path.getsize(path)
    except Exception:
        size = 0
    nrows = None
    if path.lower().endswith(".parquet") and pq is not None:
        try:
            nrows = pq.ParquetFile(path).metadata.num_rows
        except Exception:
            nrows = None
    if (nrows is not None and nrows >= cfg.auto_trigger_rows) or size >= cfg.auto_trigger_bytes:
        return cfg.auto_max_points
    return None


def _build_points_model(df_points: pd.DataFrame, coord_cols: Tuple[str, ...]) -> PointsModel:
    last_error: Optional[Exception] = None
    try:
        coords = df_points[list(coord_cols)].to_numpy(dtype=float)
        annotations = df_points.drop(columns=list(coord_cols))
        return PointsModel.parse(coords, annotations=annotations)
    except Exception as exc:
        last_error = exc
    for kwargs in [{}, {"coordinates": coord_cols}, {"coordinates": list(coord_cols)}, {"dims": coord_cols}]:
        try:
            return PointsModel.parse(df_points, **kwargs)
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Failed to build PointsModel with coord_cols={coord_cols}. Last error: {last_error}")


def load_transcripts_points(sample_dir: str, cfg: MerfishPointsConfig):
    path = find_transcripts_file(sample_dir)
    print("666666")
    print(path)
    if not path:
        return None, "no transcripts file found by patterns"
    if path.lower().endswith(".parquet"):
        if pq is None:
            return None, "pyarrow not available; cannot read parquet"
        try:
            columns = list(pq.read_schema(path).names)
        except Exception as exc:
            return None, f"cannot read parquet schema: {exc}"
    else:
        try:
            columns = list(pd.read_csv(path, nrows=0).columns)
        except Exception as exc:
            return None, f"cannot read csv header: {exc}"
    col_map = infer_transcript_columns(columns)
    if col_map["x"] is None or col_map["y"] is None or col_map["gene"] is None:
        return None, f"cannot infer required cols x/y/gene; cols head={columns[:80]}"
    cap = _auto_cap(path, cfg)
    use_cols = [col_map["x"], col_map["y"], col_map["gene"]]
    if col_map["z"] is not None:
        use_cols.append(col_map["z"])
    if cfg.keep_cell_id and col_map["cell_id"] is not None:
        use_cols.append(col_map["cell_id"])
    if cfg.keep_qv and col_map["qv"] is not None:
        use_cols.append(col_map["qv"])
    if cfg.keep_fov and col_map["fov"] is not None:
        use_cols.append(col_map["fov"])
    try:
        if path.lower().endswith(".parquet"):
            parquet_file = pq.ParquetFile(path)
            chunks = []
            loaded = 0
            for row_group in range(parquet_file.num_row_groups):
                chunk = parquet_file.read_row_group(row_group, columns=use_cols).to_pandas()
                chunks.append(chunk)
                loaded += len(chunk)
                if cap is not None and loaded >= cap:
                    break
            df = pd.concat(chunks, ignore_index=True)
            if cap is not None and len(df) > cap:
                df = df.iloc[:cap].copy()
        else:
            df = pd.read_csv(path, usecols=use_cols, nrows=cap)
    except Exception as exc:
        return None, f"failed to load transcripts: {exc}"
    df_points = pd.DataFrame()
    df_points["x"] = pd.to_numeric(df[col_map["x"]], errors="coerce")
    df_points["y"] = pd.to_numeric(df[col_map["y"]], errors="coerce")
    coord_cols: List[str] = ["x", "y"]
    if col_map["z"] is not None and col_map["z"] in df.columns:
        df_points["z"] = pd.to_numeric(df[col_map["z"]], errors="coerce")
        coord_cols.append("z")
    genes = df[col_map["gene"]]
    if genes.dtype != object:
        genes = genes.astype(str)
    df_points["gene"] = genes.astype("category")
    if cfg.keep_cell_id and col_map["cell_id"] is not None and col_map["cell_id"] in df.columns:
        df_points["cell_id"] = df[col_map["cell_id"]]
    if cfg.keep_qv and col_map["qv"] is not None and col_map["qv"] in df.columns:
        df_points["qv"] = pd.to_numeric(df[col_map["qv"]], errors="coerce")
    if cfg.keep_fov and col_map["fov"] is not None and col_map["fov"] in df.columns:
        df_points["fov"] = df[col_map["fov"]]
    df_points = df_points.dropna(subset=coord_cols).reset_index(drop=True)
    try:
        points = _build_points_model(df_points, tuple(coord_cols))
    except Exception as exc:
        return None, f"PointsModel.parse failed: {exc}"
    return points, None


def inject_transcripts_points_if_missing(sdata, sample_dir: str, cfg: Optional[MerfishPointsConfig] = None) -> List[str]:
    cfg = cfg or MerfishPointsConfig()
    warnings: List[str] = []
    try:
        existing_keys = list(getattr(sdata, "points", {}).keys())
    except Exception:
        existing_keys = []
    if cfg.points_key in existing_keys:
        return warnings
    points, error = load_transcripts_points(sample_dir, cfg)
    if points is None:
        warnings.append(f"points not loaded: {error}")
        return warnings
    try:
        sdata.points[cfg.points_key] = points
    except Exception as exc:
        warnings.append(f"failed to inject points into SpatialData: {exc}")
    return warnings
