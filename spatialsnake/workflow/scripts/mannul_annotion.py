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


def annotion(adata):
  original_clusters = adata.obs['clusters']
  print(original_clusters)
  new_categories = original_clusters.astype('string').map(cell_annotation)
  adata.obs["celltype"] = new_categories.astype('category')
  print(adata.obs["celltype"])
  return adata

if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
  adata = annotion(adata)
  sq.pl.spatial_scatter(adata, shape=None, color="celltype")
  plt.savefig(
        os.path.join(output_dir, "{args.sample_id}celltype.png"),
        dpi=300,
        bbox_inches='tight')
  plt.show()
  plt.close()
  exit()
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  for table in concatenated_sdata.tables.keys():
    table=table
    concatenated_sdata[table]=annotion(concatenated_sdata[table])
    print(concatenated_sdata,concatenated_sdata[table])





def crop0(x,crs,x1,x2,y1,y2):
    return bounding_box_query(
        x,
        min_coordinate=[x1, y1],
        max_coordinate=[x2, y2],
        axes=("x", "y"),
        target_coordinate_system=crs)
    

sample_cnt=args.sample_cnt


image_elements = list(concatenated_sdata.images.keys())
shape_elements = list(concatenated_sdata.shapes.keys())
valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)


x1,x2,y1,y2=args.coord

image_counts=len(image_elements)//sample_cnt
shape_count=len(shape_elements)//sample_cnt
shapes=[]
images=[]
systems=[]


# images = [img for img in image_elements if args.image_type in img] if image_counts=len(image_elements)//sample_cnt > 1 else image_elements
# systems = [sys for sys in valid_coord_systems if args.image_type in sys] if len(valid_coord_systems) > 1 else valid_coord_systems
# shapes = [shape for shape in shape_elements if args.shape_type in shape] if len(shape_elements)//sample_cnt > 1 else shape_elements

if image_counts>1:
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
print(args.type)
print(shape_count)
if shape_count>1 and args.type !="visium":
  for i in range(len(shape_elements)):
    if args.shape_type in shape_elements[i]:
      shapes.append(shape_elements[i])
else:
  shapes=shape_elements

print(shapes)
print(images)
print(systems)
extents=[]
for i in range(len(images)):
    sys=systems[0] if len(valid_coord_systems)<2 else systems[i]
    title=images[i].replace("_hires_image","")
    print("##########################################")
    extent = spd.get_extent(concatenated_sdata,elements=[shapes[i]],coordinate_system=sys)
    extents.append(extent)

if len(images) != len(shapes):
    print(extent)
    print("Check the spatial data to make sure that for every image there is a shape")
    print(images,shapes)
    exit()

else:
  for i in range(len(images)):
        sys=systems[0] if len(valid_coord_systems)<2 else systems[i]
        print("Plotting: "+ images[i])
        print(sys)
        print(images[i],shapes[i])
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        if args.image_slice == True:
            print(args.image_slice)
            print("sdasdasdasda")
            concatenated_sdata=crop0(concatenated_sdata,sys,x1,x2,y1,y2)
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image",coordinate_systems=sys)
        concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i],color="celltype").pl.show(ax=axes[1],coordinate_systems=sys, title=title)
        plt.savefig(
        os.path.join(dir_path, f"{images[i]}.png"),
        dpi=300,
        bbox_inches='tight')
        plt.close()
        
sc.pl.umap(concatenated_sdata[table],color="celltype",wspace=0.4)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}UMAP.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()

sc.pl.violin(concatenated_sdata[table], ["n_genes_by_counts"], groupby="celltype", rotation=90)
plt.savefig(
    os.path.join(dir_path, f"{args.sample_id}gene_enrich.png"),
      dpi=300,
      bbox_inches='tight')
plt.show()
plt.close()



if type!="slide_seq":
  concatenated_sdata[table].obs['cell_id'] = concatenated_sdata[table].obs['cell_id'].astype(str)
  concatenated_sdata[table].obs['region'] = concatenated_sdata[table].obs['region'].astype('category')
  concatenated_sdata.write(os.path.join(args.output_zarr_path),overwrite=True)
else:
  adata.write(os.path.join(args.output_zarr_path))
