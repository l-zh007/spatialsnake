import os
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import spatialdata_io
import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon

type='visium_segment'

def create_zarr_bin(path_to_outputs,zarr_name,bin_size):
    print(zarr_name)
    dataset=zarr_name.split('.')[0]
    print(dataset)
    sdata = spatialdata_io.visium_hd(path_to_outputs, dataset_id=dataset, filtered_counts_file=True, bin_size=bin_size)
    sdata.write(zarr_name, overwrite=True)
    del sdata

def creat_zarr_visium(path_to_outputs,zarr_name,h5_name):
    print(zarr_name)
    dataset=zarr_name.split('.')[0]
    sdata=spatialdata_io.visium(path_to_outputs,dataset_id=dataset,counts_file=h5_name)
    sdata.write(zarr_name, overwrite=True)
    del sdata

def create_zarr(count_matrix_path,image_path,scale_factors_path,geojson_path,sample_name):
    print(sample_name)

    # Load and Prepare Raw Data
    # Define file paths
    COUNT_MATRIX_PATH = count_matrix_path
    IMAGE_PATH = image_path
    SCALE_FACTORS_PATH = scale_factors_path
    GEOJSON_PATH = geojson_path

    # Load AnnData
    adata = sc.read_10x_h5(COUNT_MATRIX_PATH)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs.index = sample_name +"_" + adata.obs.index.astype(str)

    # Load and preprocess image data
    image_data = np.array(Image.open(IMAGE_PATH))
    if image_data.ndim == 2:
        image_data = image_data[np.newaxis, :, :] # Add channel dimension for grayscale
    elif image_data.ndim == 3:
        image_data = np.transpose(image_data, (2, 0, 1)) # (H, W, C) -> (C, H, W) for spatialdata

    # Load scale factors
    with open(SCALE_FACTORS_PATH, 'r') as f:
        scale_data = json.load(f)

    # Load GeoJSON data
    with open(GEOJSON_PATH, 'r') as f:
        geojson_data = json.load(f)

    # Define coordinate systems:
    # `downscale_to_hires`: The coordinate system where shapes are located, scaled relative to the hires resolution.

    hires_scale = scale_data['tissue_hires_scalef']

    # Transformation for shapes (from pixel to downscale_to_hires)
    shapes_transformations = {
       "downscale_to_hires": Scale(np.array([hires_scale, hires_scale]), axes=("x", "y")) # if the high-resolution microscope image is being used and Identity() transform would be performed.
    }

    # Transformation for the 'hires_tissue_image' (it's already in the 'downscale_to_hires' space visually)
    image_transformations = {
        "downscale_to_hires": Identity()
    }

    # Process Cell Segmentation (GeoJSON) and Integrate with AnnData

    # Create a mapping from adata.obs.index to geojson features
    geojson_features_map = {
        f"{sample_name}_cellid_{feature['properties']['cell_id']:09d}-1": feature
        for feature in geojson_data['features']
    }

    # Prepare data for GeoDataFrame and update adata.obs
    geometries = []
    cell_ids_ordered = []

    for obs_index_str in adata.obs.index:
        feature = geojson_features_map.get(obs_index_str)
        if feature:
            # Create shapely Polygon from coordinates
            polygon_coords = np.array(feature['geometry']['coordinates'][0])
            geometries.append(Polygon(polygon_coords))
            cell_ids_ordered.append(obs_index_str)
        else:
            geometries.append(None) # Or a suitable placeholder
            cell_ids_ordered.append(obs_index_str)

    # Remove None entries if any (or handle them upstream)
    valid_indices = [i for i, geom in enumerate(geometries) if geom is not None]
    geometries = [geometries[i] for i in valid_indices]
    cell_ids_ordered = [cell_ids_ordered[i] for i in valid_indices]


    # Create GeoDataFrame for shapes
    shapes_gdf = gpd.GeoDataFrame({
        'cell_id': cell_ids_ordered,
        'geometry': geometries
    }, index=cell_ids_ordered)
    # Update adata.obs with cluster information and spatial identifiers
    adata.obs['cell_id'] = adata.obs.index
    adata.obs['region'] = sample_name + '_cell_boundaries'
    adata.obs['region'] = adata.obs['region'].astype('category')
    adata = adata[shapes_gdf.index].copy() # Filter adata to match shapes_gdf

    # Define names for SpatialData elements
    IMAGE_KEY =  sample_name + '_hires_tissue_image'
    TABLE_KEY =  'table'
    SHAPES_KEY = sample_name + '_cell_boundaries'

    # Create SpatialData elements directly
    sdata = spd.SpatialData(
        images={
            IMAGE_KEY: Image2DModel.parse(image_data, transformations=image_transformations)
        },
        tables={
            TABLE_KEY: TableModel.parse(
                adata,
                region=SHAPES_KEY, # Link table to shapes element
                region_key='region', # Column in adata.obs indicating region name
                instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
            )
        },
        shapes={
            SHAPES_KEY: ShapesModel.parse(shapes_gdf, transformations=shapes_transformations)
        }
    )

    sdata.write(sample_name, overwrite=True)
    del sdata
    gc.collect()

def create_zarr_xenium(
                path_to_outputs,
                zarr_name
):
    print(zarr_name)
    sdata = spatialdata_io.xenium(path_to_outputs,cells_boundaries=False, nucleus_boundaries=False,cells_labels=True, nucleus_labels=False,morphology_mip=False)
    sdata.write(zarr_name,overwrite=True)
    del sdata
    


# samples = {"Non_Lesional_1":["./data/ST_21_NL","Non_Lesional_1.zarr"],
#           "Non_Lesional_2":["./data/ST_22_NL","Non_Lesional_2.zarr"]}


samples = {"Colon_Cancer_P1":["data/Cancer_P1_filtered_feature_cell_matrix.h5",
                      "data/Cancer_P1_tissue_hires_image.png",
                      "data/Cancer_P1_scalefactors_json.json",
                      "data/Cancer_P1_cell_segmentations.geojson",
                      "Colon_Cancer_P1"],
            "Colon_Cancer_P2":["data/Cancer_P2_filtered_feature_cell_matrix.h5",                       "data/Cancer_P2_tissue_hires_image.png",
                      "data/Cancer_P2_scalefactors_json.json",
                      "data/Cancer_P2_cell_segmentations.geojson",
                     "Colon_Cancer_P2"]}

# samples = {"Colon_Cancer_P1":["data/Visium_HD_Human_Colon_Cancer_P1","Colon_Cancer_P1.zarr",8],
#             "Colon_Cancer_P2":["data/Visium_HD_Human_Colon_Cancer_P2","Colon_Cancer_P2.zarr",8]}

# samples = {
#     "Kidney_Cancer":["./data/Kidney_Cancer_data","Kidney_Cancer.zarr"],
#     "Kidney_Normal":["./data/Kidney_Normal_data","Kidney_Normal.zarr"]}

samples = {
    "normal_1": ["data/normal1.h5ad","normal1.h5ad"],
    "normal_2": ["data/normal2.h5ad","normal2.h5ad"]
}


print("Saving zarr files")


if type=='visium_segment':
  for key, inputs in samples.items():
        create_zarr(count_matrix_path=inputs[0],
                image_path=inputs[1],
                scale_factors_path=inputs[2],
                geojson_path=inputs[3],
                sample_name=inputs[4])
elif type=='visium_HD':
  for key, inputs in samples.items():
        create_zarr_bin(path_to_outputs=inputs[0],
                zarr_name=inputs[2],
                bin_size=inputs[1])

  del samples, inputs, key
  gc.collect()
elif type=="visium":
  for key, inputs in samples.items():
        creat_zarr_visium(
                path_to_outputs=inputs[0],
                zarr_name=inputs[1],
                h5_name='filtered_feature_bc_matrix.h5')
elif type=="slide_seq":
  for key, inputs in samples.items():
    sdata = sc.read_h5ad(inputs[0])
    sdata.write(inputs[1])
elif type=="xenium":
  for key, inputs in samples.items():
      create_zarr_xenium(path_to_outputs=inputs[0],zarr_name=inputs[1])
  
