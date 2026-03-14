import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "8"
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import squidpy as sq
import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
import warnings
warnings.filterwarnings("ignore")
from spatialsnake.workflow.function.plot import cluster_proportion
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv

import argparse
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--image_type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--coord', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_cnt', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--markers_algorithm', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--shape_type', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--image_slice', type=str,required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()
type=args.type
dir_path=os.path.dirname(args.output_zarr_path)

def render_spatial_plots(concatenated_sdata,shapes,images,systems,extents):
    for i in range(len(images)):
        sys=systems[0] if len(valid_coord_systems)<2 else systems[i]
        title=images[i].replace("_hires_image","")
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        if args.image_slice == True:
            concatenated_sdata=crop0(concatenated_sdata,sys,x1,x2,y1,y2)
        try:
          concatenated_sdata=crop(concatenated_sdata,sys,extents[i])
        except:
          pass
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image",coordinate_systems=sys)
        if shapes[i] in ["cell_circles","cell_boundaries","nucleus_boundaries"]:
          concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i], color="clusters",table_name="table").pl.show(ax=axes[1],coordinate_systems=sys, title=title)
        else:  
          concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i],color="clusters").pl.show(ax=axes[1],coordinate_systems=sys, title=title)
        plt.savefig(
        os.path.join(dir_path, f"{images[i]}_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
        plt.close()

def crop0(x,crs,x1,x2,y1,y2):
    return bounding_box_query(
        x,
        min_coordinate=[x1, y1],
        max_coordinate=[x2, y2],
        axes=("x", "y"),
        target_coordinate_system=crs)

def crop(x,crs,bbox):
    return spd.bounding_box_query(
        x,
        min_coordinate=[bbox['x'][0], bbox['y'][0]],
        max_coordinate=[bbox['x'][1], bbox['y'][1]],
        axes=("x", "y"),
        target_coordinate_system=crs)

def slide_seq_plot(adata):
  for sample in adata.obs['region'].unique():
    adata_sample = adata[adata.obs['region'] == sample]
    sq.pl.spatial_scatter(adata_sample, color="clusters", shape=None, title=sample)
    plt.savefig(
        os.path.join(dir_path, f"{args.sample_id}_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
    plt.close()

if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
  slide_seq_plot(adata)
else:
  #args.input_dir="results/useful_results/ROI_Selection_1_cells_stats.zarr"
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
    print(adata.obs)
  spatial_attrs = adata.uns.get('spatialdata_attrs', {})
  region_name = spatial_attrs.get('region')
  instance_key = spatial_attrs.get('instance_key')
  if isinstance(region_name, list):
      region_name = region_name[0] if len(region_name) > 0 else None
  if region_name in concatenated_sdata.shapes and instance_key in adata.obs:
      adata.obs[instance_key] = adata.obs[instance_key].astype(str)
      shapes_df = concatenated_sdata.shapes[region_name]
      shapes_df.index = shapes_df.index.astype(str)
      print(concatenated_sdata.shapes[region_name])
      print("@@@@@")
      common_index = shapes_df.index.intersection(adata.obs[instance_key])
      if len(common_index) > 0:
          if len(common_index) != len(shapes_df):
              shapes_df = shapes_df.loc[common_index]
          if len(common_index) != adata.n_obs:
              adata = adata[adata.obs[instance_key].isin(common_index)].copy()
          concatenated_sdata.shapes[region_name] = shapes_df
          concatenated_sdata[table] = adata
  sample_cnt=args.sample_cnt
  image_elements = list(concatenated_sdata.images.keys())
  shape_elements = list(concatenated_sdata.shapes.keys())
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  image_counts=len(image_elements)//sample_cnt
  shape_count=len(shape_elements)//sample_cnt
  shapes=[]
  images=[]
  systems=[]
  if len(valid_coord_systems)>1:
      for i in range(len(image_elements)):
        if args.image_type in image_elements[i]:
          images.append(image_elements[i])
  else:
      images=image_elements

  if len(valid_coord_systems)>1:
      for i in range(len(valid_coord_systems)):
          if args.image_type in valid_coord_systems[i]:
            systems.append(valid_coord_systems[i])
  else:
      systems=valid_coord_systems
  print(shape_count,shape_elements)
  if shape_count>1 and type!="visium" and type != "xenium":
      for i in range(len(shape_elements)):
        if args.shape_type in shape_elements[i]:
          shapes.append(shape_elements[i])
  else:
      shapes=shape_elements
  print(f"region字段唯一值: {concatenated_sdata[table].obs['region'].unique()}")
  print(concatenated_sdata[table].uns['spatialdata_attrs'])
  print(shapes)
  print(images)
  print(systems)

  extents = []
  for i in range(len(image_elements)):
    extent =  spd.get_extent(concatenated_sdata,elements=[shape_elements[i]],coordinate_system=systems[i])
    extents.append(extent)
  render_spatial_plots(concatenated_sdata,shapes,images,systems,extents)  
    
    
    
    

cluster_proportion(adata, sample_col='region', 
                      cluster_col='clusters',
                      figsize=(12, 8),
                      palette='tab20',
                      sort_samples=None,
                      sort_clusters=None,
                      title="Clusters_proportion",
                      save_path=os.path.join(dir_path, f"Clusters_proportion.png"),
                      dpi=300)
                      
                      

markers_algorithm = args.markers_algorithm
sc.tl.rank_genes_groups(adata = adata, groupby="clusters", method=markers_algorithm)
sc.pl.rank_genes_groups_dotplot(
    adata=adata, 
    groupby="clusters", 
    standard_scale="var", 
    n_genes=5,
    show=False)
plt.savefig(os.path.join(dir_path, f"{args.sample_id}rank_genes_groups_dotplot.png"),dpi=300,bbox_inches='tight')



export_cluster_csv(
    concatenated_sdata,
    args.type,
    dir_path,
    cell_id_col="cell_id",
    info_col="clusters",
    sample_col="sample",
    sample_id=args.sample_id
)

df_marker_genes = sc.get.rank_genes_groups_df(adata = adata,group = None,pval_cutoff=0.05)
output_filename = f"marker_genes_pval.csv"
output_path = os.path.join(dir_path,output_filename)
df_marker_genes.to_csv(output_path)



for group_value, sub_df in df_marker_genes.groupby("group"):
    sub_df_filtered = sub_df.drop(columns=["group"])
    output_filename = f"cluster_{group_value}.csv"
    output_path = os.path.join(dir_path,f'{group_value}')
    os.makedirs(os.path.join(output_path), exist_ok=True)
    sub_df_filtered.to_csv(os.path.join(output_path,output_filename), index=False)
