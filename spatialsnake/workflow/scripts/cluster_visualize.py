import os
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "1"
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

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
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
parser.add_argument('--x1,x2,y1,y2', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_cnt', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--markers_algorithm', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--coord', type=int,nargs=4, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--shape_type', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--image_slice', type=str,required=False,
                   help='Path for the output zarr file')


args = parser.parse_args()

type=args.type
dir_path=os.path.dirname(args.output_zarr_path)


if type=="slide_seq":
  adata = sc.read_h5ad(args.input_dir)
  sq.pl.spatial_scatter(adata, shape=None, color="clusters")
  plt.savefig(
        os.path.join(output_dir, f"{args.sample_id}Clusters.png"),
        dpi=300,
        bbox_inches='tight')
  plt.show()
  plt.close()
  exit()
else:
  concatenated_sdata = spd.read_zarr(args.input_dir)
  print(concatenated_sdata)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]
    

# if slice:
#   def crop0(x,crs,bbox):
#     return spd.bounding_box_query(
#         x,
#         min_coordinate=[bbox['x'][0], bbox['y'][0]],
#         max_coordinate=[bbox['x'][1], bbox['y'][1]],
#         axes=("x", "y"),
#         target_coordinate_system=crs,)
# else:
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
print(type)
if shape_count>1 and type!="visium":
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
        # title=image_elements[i].replace("_hires_tissue_image","")
        title=images[i].replace("_hires_image","")
        print(f'{title}_downscaled_hires')
        print(images[i],shapes[i])
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        if args.image_slice == True:
            print(args.image_slice)
            print("sdasdasdasda")
            concatenated_sdata=crop0(concatenated_sdata,sys,x1,x2,y1,y2)
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image",coordinate_systems=sys)
        concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i],color="clusters").pl.show(ax=axes[1],coordinate_systems=sys, title=title)
        plt.savefig(
        os.path.join(dir_path, f"{images[i]}_Clusters.png"),
        dpi=300,
        bbox_inches='tight')
        plt.close()



# sq.gr.spatial_neighbors(adata, coord_type="generic")
# sq.gr.nhood_enrichment(adata, cluster_key="clusters")
# sq.pl.nhood_enrichment(adata, cluster_key="clusters", figsize=(5, 5))
# plt.savefig(
#         os.path.join(dir_path, f"{args.sample_id}Cells_Clusters.png"),
#         dpi=300,
#         bbox_inches='tight')
# plt.show()
# plt.close()


print(adata)



# t-test、wilcoxon logreg
markers_algorithm = args.markers_algorithm

sc.tl.rank_genes_groups(adata = adata, groupby="clusters", method=markers_algorithm)
sc.pl.rank_genes_groups_dotplot(
    adata=adata, 
    groupby="clusters", 
    standard_scale="var", 
    n_genes=5,
    show=False)
plt.savefig(os.path.join(dir_path, f"{args.sample_id}rank_genes_groups_dotplot.png"),dpi=300,bbox_inches='tight')

# sc.pl.stacked_violin(adata, marker_genes, groupby="clusters")
# plt.savefig(os.path.join(dir_path, f"{args.sample_id}rank_genes_groups_violin.png"),dpi=300,bbox_inches='tight')
print("Available groups in rank_genes_groups:", list(adata.uns['rank_genes_groups']['names'].dtype.names) if 'rank_genes_groups' in adata.uns else "No rank_genes_groups found")


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

# import os
# import shutil
# import argparse
# 
# def cleanup_spatialdata(input_path):
#     if os.path.exists(input_path):
#         print(f"Cleaning up spatialdata: {input_path}")
#         if os.path.isdir(input_path):
#             shutil.rmtree(input_path)
#         else:
#             os.remove(input_path)
#     else:
#         print(f"Spatialdata path does not exist: {input_path}")
# cleanup_spatialdata(args.input_dir)
# 
# 
# if type!="slide_seq":
#   concatenated_sdata[table]=adata
#   concatenated_sdata.write(os.path.join(args.input_dir),overwrite=True)
# else:
#   adata.write(os.path.join(args.input_dir))












# print(image_elements,shape_elements)
# print(concatenated_sdata.shapes)
# print(concatenated_sdata)
# print(concatenated_sdata["segmentation_counts"])
# print("##########################################")
# print(concatenated_sdata)
# print(f"region字段唯一值: {concatenated_sdata['segmentation_counts'].obs['region'].unique()}")
# print(concatenated_sdata['segmentation_counts'].obs)
# print(concatenated_sdata['segmentation_counts'].uns['spatialdata_attrs'])





# sq.gr.spatial_neighbors(adata, coord_type="generic")
# sq.gr.nhood_enrichment(adata, cluster_key="clusters")
# sq.pl.nhood_enrichment(adata, cluster_key="clusters", figsize=(5, 5))
# plt.savefig(
#         os.path.join(output_dir, f"{key}Cells_Clusters.png"),
#         dpi=300,
#         bbox_inches='tight')
# plt.show()
# plt.close()








