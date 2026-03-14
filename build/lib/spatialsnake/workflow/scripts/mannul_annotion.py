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
import argparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
from spatialsnake.workflow.function.plot import cluster_proportion
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
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
parser.add_argument('--sample_cnt', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--coord', type=int,nargs=4, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--anno_data', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--image_slice', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--shape_type', type=str,required=False,
                   help='Path for the output zarr file')

args = parser.parse_args()
dir_path=os.path.dirname(args.output_zarr_path)
if args.anno_data:
  print(args.anno_data,args.sample_id)
  cell_annotation = json.loads(args.anno_data)
  cell_annotation = cell_annotation[args.sample_id]
  print(cell_annotation,args.sample_id)
  print("####################")
# cell_annotation = {
#     '0': 'Fibroblast',
#     '1': 'Tumor',
#     '2': 'Plasma cell',
#     '3': 'Smooth muscle',
#     '4': 'Tumor',
#     '5': 'B cell and T cell',
#     '6': 'Enterocyte',
#     '7': 'Plasma cell',
#     '8': 'Goblet',
#     '9': 'Enterocyte'
# }
def render_spatial_plots(concatenated_sdata,shapes,images,systems):
    print(images,shapes,systems,"@@@@@@@@@@@@@")
    for i in range(len(images)):
        sys=systems[0] if len(valid_coord_systems)<2 else systems[i]
        title=images[i].replace("_hires_image","")
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        if args.image_slice == True:
            concatenated_sdata=crop0(concatenated_sdata,sys,x1,x2,y1,y2)
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image",coordinate_systems=sys)
        concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i],color="celltype").pl.show(ax=axes[1],coordinate_systems=sys, title=title)
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

def slide_seq_plot(adata):
  for sample in adata.obs['region'].unique():
    adata_sample = adata[adata.obs['region'] == sample]
    sq.pl.spatial_scatter(adata_sample, color="celltype", shape=None, title=sample)
    plt.savefig(
        os.path.join(dir_path, f"{args.sample_id}_celltype.png"),
        dpi=300,
        bbox_inches='tight')
    plt.close()

def annotion(adata,cell_annotation):
  original_clusters = adata.obs['clusters']
  print(original_clusters)
  new_categories = original_clusters.astype('string').map(cell_annotation)
  adata.obs["celltype"] = new_categories.astype('category')
  print(adata.obs["celltype"])
  return adata

if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
  adata = annotion(adata,cell_annotation)
  slide_seq_plot(adata)
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
  adata = annotion(adata,cell_annotation)
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
  #if shape_count>1 and type!="visium":
      #for i in range(len(shape_elements)):
        #if args.shape_type in shape_elements[i]:
          #shapes.append(shape_elements[i])
      #print(shapes,"66666")
  shapes=shape_elements
  render_spatial_plots(concatenated_sdata,shapes,images,systems) 


cluster_proportion(adata, sample_col='region', 
                      cluster_col='celltype',
                      figsize=(12, 8),
                      palette='tab20',
                      sort_samples=None,
                      sort_clusters=None,
                      title="celltype_proportion",
                      save_path=os.path.join(dir_path, f"celltype_proportion.png"),
                      dpi=300)
        
sc.pl.umap(adata,color="celltype",wspace=0.4)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}UMAP.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()

sc.pl.violin(adata, ["n_genes_by_counts"], groupby="celltype", rotation=90)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}_gene_enrich.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()

export_cluster_csv(
  adata if args.type=="slide_seq" else concatenated_sdata,
  args.type,
  dir_path,
  cell_id_col="cell_id",
  info_col="celltype",
  sample_col="region",
  sample_id=args.sample_id
)



if type!="slide_seq":
  concatenated_sdata[table].obs['cell_id'] = concatenated_sdata[table].obs['cell_id'].astype(str)
  concatenated_sdata[table].obs['region'] = concatenated_sdata[table].obs['region'].astype('category')
  concatenated_sdata.write(os.path.join(args.output_zarr_path),overwrite=True)
else:
  adata.write(os.path.join(args.output_zarr_path))
