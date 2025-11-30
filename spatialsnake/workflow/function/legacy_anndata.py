import warnings
from typing import TYPE_CHECKING

import numpy as np
from spatialdata import (
    SpatialData,
    get_centroids,
    get_extent,
    join_spatialelement_table,
    to_circles,
)
from spatialdata._core.operations._utils import transform_to_data_extent
from spatialdata.models import Image2DModel, ShapesModel, TableModel, get_table_keys
from spatialdata.transformations import Identity, Scale

from anndata import AnnData



def transform_to_zarr(adata: AnnData,sample_id: str) -> SpatialData:
    """Convert (legacy) spatial AnnData object to SpatialData object.
    
     revise from  https://github.com/scverse/spatialdata-io/blob/main/src/spatialdata_io/converters/legacy_anndata.py
     Add revice some data to match spatialsnake pipeline
    """
    SPATIAL = "spatial"

    SCALEFACTORS = "scalefactors"
    TISSUE_HIRES_SCALEF = "tissue_hires_scalef"
    TISSUE_LOWRES_SCALEF = "tissue_lowres_scalef"
    SPOT_DIAMETER_FULLRES = "spot_diameter_fullres"

    IMAGES = "images"
    HIRES = "hires"
    LOWRES = "lowres"

    # SpatialData keys
    REGION = sample_id
    REGION_KEY = "region"
    INSTANCE_KEY = "cell_id"
    SPOT_DIAMETER_FULLRES_DEFAULT = 10

    images = {}
    shapes = {}
    spot_diameter_fullres_list = []
    shapes_transformations = {}

    if SPATIAL in adata.uns:
        dataset_ids = list(adata.uns[SPATIAL].keys())
        for dataset_id in dataset_ids:
            # read the image data and the scale factors for the shapes
            keys = set(adata.uns[SPATIAL][dataset_id].keys())
            tissue_hires_scalef = None
            tissue_lowres_scalef = None
            hires = None
            lowres = None
            if SCALEFACTORS in keys:
                scalefactors = adata.uns[SPATIAL][dataset_id][SCALEFACTORS]
                if TISSUE_HIRES_SCALEF in scalefactors:
                    tissue_hires_scalef = scalefactors[TISSUE_HIRES_SCALEF]
                if TISSUE_LOWRES_SCALEF in scalefactors:
                    tissue_lowres_scalef = scalefactors[TISSUE_LOWRES_SCALEF]
                if SPOT_DIAMETER_FULLRES in scalefactors:
                    spot_diameter_fullres_list.append(scalefactors[SPOT_DIAMETER_FULLRES])
            if IMAGES in keys:
                image_data = adata.uns[SPATIAL][dataset_id][IMAGES]
                if HIRES in image_data:
                    hires = image_data[HIRES]
                if LOWRES in image_data:
                    lowres = image_data[LOWRES]

            # construct the spatialdata elements
            if hires is not None:
                if 'hires_image' in dataset_id:
                  image_index=f"{dataset_id}"
                else:
                  image_index=f"{dataset_id}_hires_image"
                # prepare the hires image
                assert tissue_hires_scalef is not None, (
                    "tissue_hires_scalef is required when an the hires image is present"
                )
                hires_image = Image2DModel.parse(
                    hires, dims=("y", "x", "c"), transformations={f"{image_index}": Identity()}
                )
                images[image_index] = hires_image

                # prepare the transformation to the hires image for the shapes
                scale_hires = Scale([tissue_hires_scalef, tissue_hires_scalef], axes=("x", "y"))
                shapes_transformations[image_index] = scale_hires
            if lowres is not None:
                if 'hires_image' in dataset_id:
                  image_index=f"{dataset_id}"
                else:
                  image_index=f"{dataset_id}_lowres_image"
                assert tissue_lowres_scalef is not None, (
                    "tissue_lowres_scalef is required when an the lowres image is present"
                )
                lowres_image = Image2DModel.parse(
                    lowres, dims=("y", "x", "c"), transformations={f"{image_index}": Identity()}
                )
                images[image_index] = lowres_image

                # prepare the transformation to the lowres image for the shapes
                scale_lowres = Scale([tissue_lowres_scalef, tissue_lowres_scalef], axes=("x", "y"))
                shapes_transformations[image_index] = scale_lowres

    # validate the spot_diameter_fullres value
    if len(spot_diameter_fullres_list) > 0:
        d = np.array(spot_diameter_fullres_list)
        if not np.allclose(d, d[0]):
            warnings.warn(
                "spot_diameter_fullres is not constant across datasets. Using the average value.",
                UserWarning,
                stacklevel=2,
            )
            spot_diameter_fullres = d.mean()
        else:
            spot_diameter_fullres = d[0]
    else:
        warnings.warn(
            f"spot_diameter_fullres is not present. Using {SPOT_DIAMETER_FULLRES_DEFAULT} as default value.",
            UserWarning,
            stacklevel=2,
        )
        spot_diameter_fullres = SPOT_DIAMETER_FULLRES_DEFAULT

    # parse and prepare the shapes
    if SPATIAL in adata.obsm:
        xy = adata.obsm[SPATIAL]
        radius = spot_diameter_fullres / 2
        shapes[REGION] = ShapesModel.parse(xy, geometry=0, radius=radius, transformations=shapes_transformations)
        # link the shapes to the table
        adata.obs['cell_id'] = adata.obs.index
        new_table = adata.copy()
        if TableModel.ATTRS_KEY in new_table.uns:
            del new_table.uns[TableModel.ATTRS_KEY]
        new_table.obs[REGION_KEY] = REGION
        new_table.obs[REGION_KEY] = new_table.obs[REGION_KEY].astype("category")
        new_table.obs[INSTANCE_KEY] = shapes[REGION].index.values
        new_table = TableModel.parse(new_table, region=REGION, region_key=REGION_KEY, instance_key=INSTANCE_KEY)
    else:
        new_table = adata.copy()
    return SpatialData(tables={"table": new_table}, images=images, shapes=shapes)
