import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import sys
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
import argparse
import matplotlib as mpl
import cell2location
from matplotlib import rcParams

from typing import Literal
from spatialdata.datasets import blobs_annotating_element
from spatialdata.transformations import Affine, set_transformation
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
from cell2location.utils.filtering import filter_genes
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", required=True)
parser.add_argument("--sample_id", required=True)
parser.add_argument("--output_zarr_path", required=True)
parser.add_argument("--type",  required=True)
parser.add_argument("--image_type", required=True)
parser.add_argument("--image_slice", required=False)
parser.add_argument("--sample_cnt",type=int, required=True)
parser.add_argument('--coord', type=int,nargs=4, required=False,
                   help='Path for the output zarr file')
parser.add_argument("--input_st", required=False)
parser.add_argument("--shape_type", required=False)
args = parser.parse_args()




print(args.image_slice)
print(args.output_zarr_path)
output_dir=os.path.dirname(args.output_zarr_path)
os.makedirs(os.path.join(output_dir,"figure"), exist_ok=True)

concatenated_sdata = spd.read_zarr(args.input_dir)
for table in concatenated_sdata.tables.keys():
    table=table
    adata_vis = concatenated_sdata[table]
if args.input_st=="":
  flag=False
else:
  flag=True
  input_st=os.path.dirname(args.input_st)
  adata = sc.read_visium(input_st, count_file='filtered_feature_bc_matrix.h5', load_images=True)
# print(adata)
# adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
# adata_vis.var['SYMBOL'] = adata_vis.var_names
# adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# adata_file = "./sp.h5ad"
# adata_vis = sc.read_h5ad(adata_file)


# mod = cell2location.models.Cell2location.load(args.input_dir, adata_vis)
# concatenated_sdata[table] = mod.export_posterior(
#     concatenated_sdata[table], sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
# )
# 
# 
# 
# mod.plot_QC()
# plt.savefig('04-mod.plot_QC.png')

# print(adata_vis.obsm['q05_cell_abundance_w_sf'])
# adata_vis.obsm['q05_cell_abundance_w_sf']
# adata_vis.obsm['q95_cell_abundance_w_sf']


concatenated_sdata[table].obs[concatenated_sdata[table].uns['mod']['factor_names']] = concatenated_sdata[table].obsm['q05_cell_abundance_w_sf']

print(concatenated_sdata)


# select one slide
# from cell2location.utils import select_slide
# slide = select_slide(concatenated_sdata[table], 'Human_Lymph_Node')
# print("可用的 img_key 列表：", slide.uns['spatial'].keys())



print(concatenated_sdata[table])
print(concatenated_sdata[table].uns["mod"]["factor_names"])




new_columns = concatenated_sdata[table].obs.columns.str.replace("+", "plus", regex=False)
concatenated_sdata[table].obs.columns = new_columns
print("\n修改后的 obs 列名：")
print(concatenated_sdata[table].obs.columns.tolist())

new_col=[]
for i in concatenated_sdata[table].uns["mod"]["factor_names"]:
  new_col.append(i.replace("+", "plus"))




existing_cell_cols = [col for col in new_col if col in concatenated_sdata[table].obs.columns]
if len(existing_cell_cols) != len(new_col):
    missing_cols = set(new_col) - set(existing_cell_cols)
    print(f"警告：以下细胞类型列在 obs 中不存在，已自动过滤：{missing_cols}")

abundance_matrix = concatenated_sdata[table].obs[existing_cell_cols]

concatenated_sdata[table].obs["cellLoca_type"] = abundance_matrix.idxmax(axis=1)

print(concatenated_sdata[table].obs["cellLoca_type"].head())








x1,x2,y1,y2=args.coord
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
        print(images[i],shapes[i])
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        if args.image_slice == True:
            print(args.image_slice)
            concatenated_sdata=crop0(concatenated_sdata,sys,x1,x2,y1,y2)
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image",coordinate_systems=sys)
        concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i], color='cellLoca_type').pl.show(ax=axes[1],coordinate_systems=sys)
        plt.savefig(f'{output_dir}/figure/Most_rank_celltype.png',dpi=300,bbox_inches='tight')

concatenated_sdata.write(output_dir +f"/{args.sample_id}.zarr")




############# each type 可视化
from cell2location.plt import plot_spatial
if flag:
  concatenated_sdata[table].uns["spatial"]=adata.uns["spatial"]
  sc.pl.spatial(
    concatenated_sdata[table],
    color=existing_cell_cols,
    show=False,img_key='hires',
    save=(f"{output_dir}/figure/each_celltype.png"))
  clust_labels = new_col[:6] if len(new_col)>7 else new_col
  clust_col = ['' + str(i) for i in clust_labels] 
  with mpl.rc_context({'figure.figsize': (15, 15)}):
      fig = plot_spatial(
        adata=concatenated_sdata[table],
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right')
  plt.savefig(f"{output_dir}/figure/enrich_celltype.png")



# slide = select_slide(adata_vis, 'Human_Lymph_Node')
