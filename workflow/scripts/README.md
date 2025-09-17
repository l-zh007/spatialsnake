import os
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

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
from sklearn.cluster import KMeans
import squidpy as sq
from spatialdata import bounding_box_query
# 
# def crop0(x,crs,bbox):
#     return spd.bounding_box_query(
#         x,
#         min_coordinate=[bbox['x'][0], bbox['y'][0]],
#         max_coordinate=[bbox['x'][1], bbox['y'][1]],
#         axes=("x", "y"),
#         target_coordinate_system=crs,
#     )

def crop0(x,crs,x1,x2,y1,y2):
    return bounding_box_query(
        x,
        min_coordinate=[16_000, 6000],
        max_coordinate=[22_000, 8500],
        axes=("x", "y"),
        target_coordinate_system=crs,
    )
   
samples = {"Non_Lesional_1":["./data/ST_21_NL","Non_Lesional_1.zarr"],
          "Non_Lesional_2":["./data/ST_22_NL","Non_Lesional_2.zarr"]}


# samples = {"Colon_Cancer_P1":["data/Visium_HD_Human_Colon_Cancer_P1",8,"Colon_Cancer_P1.zarr"],
#             "Colon_Cancer_P2":["data/Visium_HD_Human_Colon_Cancer_P2",8,"Colon_Cancer_P2.zarr"]}
slice=False
type="xenium"
image="hires"
shape="008um"

# samples = {
#     "normal_1": ["data/normal1.h5ad","normal_1.h5ad"],
#     "normal_2": ["data/normal2.h5ad","normal_2.h5ad"]
# }

# samples = {
#     "Kidney_Cancer":["./data/Kidney_Cancer_data","Kidney_Cancer.zarr"],
#     "Kidney_Normal":["./data/Kidney_Normal_data","Kidney_Normal.zarr"]}
pre_dir="clustering_results"
output_dir = "clustering_results"
os.makedirs(output_dir, exist_ok=True)


for key, inputs in samples.items():
    if type=="slide_seq":
      print(os.path.join(pre_dir,inputs[-1]))
      adata = sc.read_h5ad(os.path.join(pre_dir,inputs[-1]))
      sq.pl.spatial_scatter(adata, shape=None, color="clusters")
      plt.savefig(
        os.path.join(output_dir, f"{key}Clusters.png"),
        dpi=300,
        bbox_inches='tight')
      plt.show()
      plt.close()
    else:
      concatenated_sdata = spd.read_zarr(os.path.join(pre_dir,inputs[-1]))
      for table in concatenated_sdata.tables.keys():
        print(table)
        adata = concatenated_sdata[table]
        
      image_elements = list(concatenated_sdata.images.keys())
      shape_elements = list(concatenated_sdata.shapes.keys())
      valid_coord_systems = concatenated_sdata.coordinate_systems
      
      print(image_elements,shape_elements)
      print(valid_coord_systems)
      
      
      
      if len(shape_elements)>1:
         for j in range(len(shape_elements)):
           if shape in shape_elements[j]:
             shapes=shape_elements[j]
      else:
          shapes=shape_elements[0]
             
        
      if len(image_elements)>1:
         for j in range(len(image_elements)): 
            if image in image_elements[j]:
              images=image_elements[j]
      else:
         images=image_elements[0]
      if len(valid_coord_systems)>1:
        for i in range(len(valid_coord_systems)):
          if image in valid_coord_systems[i]:
            system=valid_coord_systems[i]
          else:
            system=valid_coord_systems[0]
      else:
        system=valid_coord_systems[0]
      print(image,shape)
      print(system)
     
      if slice:
        concatenated_sdata=crop0(concatenated_sdata,system,x1,x2,y1,y2)
      axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
      concatenated_sdata.pl.render_images(images).pl.show(ax=axes[0], title="image",coordinate_systems=system)
      concatenated_sdata.pl.render_images(images).pl.render_shapes(shapes,color="clusters").pl.show(ax=axes[1],title="shape",coordinate_systems=system)
      plt.savefig(
      os.path.join(output_dir, f"{key}_Clusters.png"),
      dpi=300,
      bbox_inches='tight')
      plt.close()




      sq.gr.co_occurrence(adata, cluster_key="clusters")
      sq.pl.co_occurrence(
        adata,
        cluster_key="clusters",
        clusters="0",
        figsize=(8, 4),
        )

      plt.savefig(
          os.path.join(output_dir, f"{key}co_occur.png"),
          dpi=300,
          bbox_inches='tight')
      plt.show()
      plt.close()
      
      
      sq.gr.ligrec(
      adata,
      n_perms=100,
      cluster_key="clusters",
      use_raw=False)
      
      
      
      
      sq.pl.ligrec(
        adata,
        cluster_key="clusters",
        source_groups="0",
        target_groups=["1", "2"],
        means_range=(3, np.inf),
        alpha=1e-4,
        swap_axes=True)


      plt.savefig(
          os.path.join(output_dir, f"{key}co_occur2.png"),
          dpi=300,
          bbox_inches='tight')
      plt.show()
      plt.close()
      
      
      
      
      
      
      
      

    sq.gr.spatial_neighbors(adata, coord_type="generic")
    sq.gr.nhood_enrichment(adata, cluster_key="clusters")
    sq.pl.nhood_enrichment(adata, cluster_key="clusters", figsize=(5, 5))
    plt.savefig(
        os.path.join(output_dir, f"{key}Cells_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
    plt.show()
    plt.close()

    mode = "L"
    sq.gr.ripley(adata, cluster_key="clusters", mode=mode, max_dist=500)
    sq.pl.ripley(adata, cluster_key="clusters", mode=mode)
    plt.savefig(
        os.path.join(output_dir, f"{key}Across_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
    plt.show()
    plt.close()
    
    
    
    
    
