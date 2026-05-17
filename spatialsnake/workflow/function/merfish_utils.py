from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from spatialdata.transformations import Affine, Identity, set_transformation


def normalize_optional_arg(value):
    if value is None:
        return None
    text = str(value).strip()
    if text.lower() in {"", "false", "none", "null"}:
        return None
    return text


def pick_preferred_key(keys: List[str], preferred_keywords: List[str]) -> Optional[str]:
    if not keys:
        return None
    lowered = [(key, key.lower()) for key in keys]
    for keyword in preferred_keywords:
        for original, lower in lowered:
            if keyword in lower:
                return original
    return keys[0]


def find_merfish_transform_csv(sample_dir: str) -> Optional[str]:
    patterns = [
        "images/micron_to_mosaic_pixel_transform.csv",
        "**/micron_to_mosaic_pixel_transform.csv",
        "**/*mosaic*transform*.csv",
        "**/*transform*.csv",
    ]
    candidates: List[str] = []
    for pattern in patterns:
        candidates.extend(glob.glob(os.path.join(sample_dir, pattern), recursive=True))
    return _pick_largest_file(_dedup(candidates))


def _dedup(values: List[str]) -> List[str]:
    return list(dict.fromkeys(values))


def _pick_largest_file(paths: List[str]) -> Optional[str]:
    files = [p for p in paths if os.path.isfile(p)]
    if not files:
        return None
    files.sort(key=lambda p: os.path.getsize(p), reverse=True)
    return files[0]


def load_transform_matrix(transform_csv: str) -> np.ndarray:
    transform_path = Path(transform_csv)
    for delimiter in (",", None):
        try:
            matrix = np.loadtxt(transform_path, dtype=float, delimiter=delimiter)
            if matrix.shape == (3, 3):
                return matrix
        except Exception:
            continue
    raise ValueError(f"Expected 3x3 transform matrix at {transform_csv}")


def get_scale0_coords(image) -> Tuple[float, float, float, float]:
    if hasattr(image, "__getitem__") and "scale0" in image:
        node = image["scale0"]
        dataset = node.ds if hasattr(node, "ds") else node
    else:
        dataset = image
    x = np.asarray(dataset.coords["x"].values, dtype=float)
    y = np.asarray(dataset.coords["y"].values, dtype=float)
    return float(x.min()), float(y.min()), float(x.max()), float(y.max())


def transform_points(matrix: np.ndarray, points: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    ones = np.ones((pts.shape[0], 1), dtype=float)
    homo = np.concatenate([pts, ones], axis=1)
    out = (matrix @ homo.T).T
    return out[:, :2]


def transform_bounds(matrix: np.ndarray, bounds) -> np.ndarray:
    xmin, ymin, xmax, ymax = map(float, bounds)
    corners = np.array(
        [[xmin, ymin], [xmax, ymin], [xmin, ymax], [xmax, ymax]],
        dtype=float,
    )
    transformed = transform_points(matrix, corners)
    return np.array(
        [
            transformed[:, 0].min(),
            transformed[:, 1].min(),
            transformed[:, 0].max(),
            transformed[:, 1].max(),
        ],
        dtype=float,
    )


def box_to_box_affine(src_bounds, dst_bounds) -> np.ndarray:
    sx0, sy0, sx1, sy1 = map(float, src_bounds)
    dx0, dy0, dx1, dy1 = map(float, dst_bounds)
    sw = sx1 - sx0
    sh = sy1 - sy0
    dw = dx1 - dx0
    dh = dy1 - dy0
    if sw == 0 or sh == 0:
        raise ValueError("Source bounds have zero width or height")
    ax = dw / sw
    ay = dh / sh
    tx = dx0 - ax * sx0
    ty = dy0 - ay * sy0
    return np.array(
        [[ax, 0.0, tx], [0.0, ay, ty], [0.0, 0.0, 1.0]],
        dtype=float,
    )


def shape_bounds(sdata, shape_key: str) -> np.ndarray:
    gdf = sdata.shapes[shape_key]
    if hasattr(gdf, "compute"):
        gdf = gdf.compute()
    return np.asarray(gdf.geometry.total_bounds, dtype=float)


def align_merfish_image_to_shape_space(
    sdata,
    transform_csv: str,
    image_key: Optional[str] = None,
    shape_key: Optional[str] = None,
    point_keys: Optional[List[str]] = None,
    coordinate_system: str = "global",
    register_points: bool = True,
):
    image_keys = list(getattr(sdata, "images", {}).keys())
    shape_keys = list(getattr(sdata, "shapes", {}).keys())
    if not image_keys or not shape_keys:
        raise RuntimeError("MERFISH alignment requires both image and shape elements")

    image_key = normalize_optional_arg(image_key) or pick_preferred_key(image_keys, ["z3", "mosaic", "image", "z2", "z1"])
    shape_key = normalize_optional_arg(shape_key) or pick_preferred_key(shape_keys, ["polygon", "cell", "boundar", "shape"])
    if register_points:
        point_keys = point_keys or list(getattr(sdata, "points", {}).keys())
    else:
        point_keys = []

    um_to_px = load_transform_matrix(transform_csv)
    px_to_um = np.linalg.inv(um_to_px)

    image_px_bounds = np.array(get_scale0_coords(sdata.images[image_key]), dtype=float)
    image_um_bounds = transform_bounds(px_to_um, image_px_bounds)
    local_shape_bounds = shape_bounds(sdata, shape_key)

    um_to_local = box_to_box_affine(image_um_bounds, local_shape_bounds)
    px_to_local = um_to_local @ px_to_um

    image_transform = Affine(px_to_local, input_axes=("x", "y"), output_axes=("x", "y"))
    points_transform = Affine(um_to_local, input_axes=("x", "y"), output_axes=("x", "y"))

    set_transformation(sdata.images[image_key], image_transform, to_coordinate_system=coordinate_system)
    set_transformation(sdata.shapes[shape_key], Identity(), to_coordinate_system=coordinate_system)

    for point_key in point_keys:
        if point_key in getattr(sdata, "points", {}):
            set_transformation(sdata.points[point_key], points_transform, to_coordinate_system=coordinate_system)

    return sdata
