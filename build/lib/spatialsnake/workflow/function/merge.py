import os
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt
from spatialdata.transformations import Identity, Scale
import argparse
# parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
# parser.add_argument('--INPUT', nargs='+', required=True, 
#                    help='Path to the raw data directory')
# parser.add_argument('--output_dir', type=str, required=True,
#                    help='Path for the output zarr file')
# parser.add_argument('--merge_by', type=str, required=True,
#                    help='Path for the output zarr file')
# parser.add_argument('--reordering', type=str, required=True,
#                    help='Path for the output zarr file')
# parser.add_argument('--re_sample', type=str, required=True,
#                    help='Path for the output zarr file')
# args = parser.parse_args()

def add_sample_barcode(adata,file_path):
  if 'sample' not in adata.obs.columns:
      file_name = os.path.basename(file_path)
      sample_barcode = os.path.splitext(file_name)[0]
      adata.obs['sample'] = sample_barcode
  return adata


def reorder_by_cluster(table, cluster_key, site_counter):
    if cluster_key not in table.obs.columns:
        return False, site_counter, None
    raw_order = table.obs[cluster_key].unique()
    mapping = {val: idx + site_counter for idx, val in enumerate(raw_order)}
    table.obs['clusters'] = table.obs[cluster_key].map(mapping)
    return True, site_counter + len(raw_order), table


def merge_by_samples(INPUT_list, re_sample):
    sdatas = []
    for file_path in INPUT_list:
        sdata = spd.read_zarr(file_path)
        TABLE_KEY = list(sdata.tables.keys())[0]
        table = sdata.tables[TABLE_KEY]
        if re_sample:
            table = add_sample_barcode(table, file_path)
        table.obs['cell_id'] = table.obs.index
        table.uns['spatialdata_attrs']['instance_key'] = 'cell_id'
        sdata.tables[TABLE_KEY] = table
        sdatas.append(sdata)
    return sdatas


def merge_by_clusters(INPUT_list, reordering, cluster_key):
    sdatas = []
    site_counter = 0
    for file_path in INPUT_list:
        sdata = spd.read_zarr(file_path)
        TABLE_KEY = list(sdata.tables.keys())[0]
        table = sdata.tables[TABLE_KEY]
        print(table)
        if reordering and cluster_key in table.obs.columns:
            success, site_counter, updated_table = reorder_by_cluster(
                table, cluster_key, site_counter
            )
            if not success:
                print(f"'{cluster_key}' not in the table")
                sys.exit(1)
            table = updated_table
        table.obs['cell_id'] = table.obs.index
        table.uns['spatialdata_attrs']['instance_key'] = 'cell_id'
        sdata.tables[TABLE_KEY] = table
        sdatas.append(sdata)
    return sdatas





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge')
    parser.add_argument('--INPUT', nargs='+', required=True, 
                       help='INPUT zarr')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='output_dir')
    parser.add_argument('--merge_by', type=str, 
                       required=True)
    parser.add_argument('--reordering', type=str, default='False',
                       help='重排聚类标签')
    parser.add_argument('--re_sample', type=str, default='False',
                       help='样本标识')
    parser.add_argument('--cluster_key', type=str, default='leiden',
                       help='聚类列名')
    args = parser.parse_args()
    print("starting .......................")
    print(args.INPUT)
    reordering = args.reordering.lower() == 'true'
    re_sample = args.re_sample.lower() == 'true'
    if args.merge_by == "sample":
        sdatas = merge_by_samples(args.INPUT, re_sample)
    else:
        sdatas = merge_by_clusters(args.INPUT, reordering, args.cluster_key)
    concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
    concatenated_sdata.write(os.path.join(args.output_dir,"concatenated_sdata.zarr"), overwrite=True)
    del concatenated_sdata, sdatas
    gc.collect()


# reordering = args.reordering=='True'
# re_sample = args.reordering=='True'
# if args.merge_by=="sample":
#   sdatas=merge_by_samples(args.INPUT,re_sample)
# else:
#   sdatas=merge_by_clusters(args.INPUT,reordering,barcode)
# concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
# concatenated_sdata.write(args.output_zarr_path, overwrite=True)
# del concatenated_sdata, sdatas
# gc.collect()









# if type=='visium_segment':
#   for i in range(len(args.INPUT)):
#     sdata=spd.read_zarr(args.input_path[i])
#     for table in sdata.tables.values():
#       table.var_names_make_unique()
#       # table.obs["sample"]=sample[i]
#     sdatas.append(sdata)
#     print(sdata)
#     del sdata,table
# elif type=='visium_HD':
#   for i in range(len(args.INPUT)):
#     sdata=spd.read_zarr(args.input_path[i])
#     TABLE_KEY = 'segmentation_counts'
#     for shapes in sdata.shapes.keys():
#       shape_key=shapes
#     for table in sdata.tables.values():
#         table.obs['cell_id'] = table.obs.index
#         table.obs["sample"]=sample[i]
#         table.obs["group"]=group[i]
#         table.obs['region']=shape_key
#         print(table.obs['region'])
#         print(table.uns['spatialdata_attrs'])
#         sdata.shapes[shape_key].index=table.obs.index
#     del table.uns['spatialdata_attrs']
#     sdata.tables={
#             TABLE_KEY: TableModel.parse(
#                 table,
#                 region=shape_key, # Link table to shapes element
#                 region_key='region', # Column in adata.obs indicating region name
#                 instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
#             )
#         }
#     sdatas.append(sdata)
#     print(sdata)
# elif type=="visium":
#   for i in range(len(args.INPUT)):
#       sdata=spd.read_zarr(args.input_path[i])
#       print(sdata)
#       SHAPES_KEY = sample[i]
#       TABLE_KEY = 'segmentation_counts'
#       for table in sdata.tables.values():
#           table.obs["sample"] = sample[i]
#           table.obs["group"]=group[i]
#           table.obs['cell_id'] = table.obs.index
#           sdata.shapes[sample[i]].index=table.obs['cell_id']
#           print(sdata.shapes)
#       del table.uns['spatialdata_attrs']
#       sdata.tables={
#               TABLE_KEY: TableModel.parse(
#                   table,
#                   region=SHAPES_KEY, # Link table to shapes element
#                   region_key='region', # Column in adata.obs indicating region name
#                   instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
#               )
#           }
#       print(sdata.shapes)
#       sdatas.append(sdata)
# elif type=="slide_seq":
#   for i in range(len(args.INPUT)):
#     sdata = sc.read_h5ad(args.input_path[i])
#     sdata.obs['sample'] = sample[i]
#     sdata.obs["group"]=group[i]
#     sdata.var_names_make_unique()
#     sdata.obs_names_make_unique()
#     print(sdata)
#     sdatas.append(sdata)
#   adata = anndata.concat(sdatas, join='inner', index_unique=None)
#   adata.obs_names_make_unique
#   adata.write("./concatenated_sdata")
#   exit()
# elif type=="xenium":
#   for i in range(len(args.INPUT)):
#     sdata=spd.read_zarr(args.input_path[i])
#     new_images = {}
#     for img_name in sdata.images.keys():
#         new_name = f"{sample[i]}_{img_name}"
#         new_images[new_name] = sdata.images[img_name]
#     sdata.images = new_images
#     new_images = {}
#     for shapes_name in sdata.shapes.keys():
#         new_name = f"{sample[i]}_{shapes_name}"
#         new_images[new_name] = sdata.shapes[shapes_name]
#     sdata.shapes = new_images
#     new_images = {}
#     for labels_name in sdata.labels.keys():
#         new_name = f"{sample[i]}_{labels_name}"
#         new_images[new_name] = sdata.labels[labels_name]
#     sdata.labels = new_images
# 
#     new_images = {}
#     for points_name in sdata.points.keys():
#         new_name = f"{sample[i]}_{points_name}"
#         new_images[new_name] = sdata.points[points_name]
#     sdata.points = new_images
# 
#     new_images = {}
#     SHAPES_KEY = f"{sample[i]}_{shapes_name}"
#     TABLE_KEY = 'segmentation_counts'
#     for table in sdata.tables.values():
#         table.var_names_make_unique()
#         table.obs["sample"] = sample[i]
#         table.obs["group"]=group[i]
#         # table.obs['cell_id'] = table.obs['cell_id'].astype(str)
#         table.obs['region'] = SHAPES_KEY
#     print(table.obs['region'])
#     del table.uns['spatialdata_attrs']
#     sdata.tables={
#             TABLE_KEY: TableModel.parse(
#                 table,
#                 region=SHAPES_KEY,
#                 region_key='region',
#                 instance_key='cell_id'
#             )
#         }
#     sdatas.append(sdata)










