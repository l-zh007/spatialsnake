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
import anndata
import argparse
import seaborn as sns
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_path', nargs='+', required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--group', nargs='+', required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', nargs='+', required=False,
                   help='Path for the output zarr file')
args = parser.parse_args()

type=args.type
group=args.group
sample=args.sample_id
def QC_plot(type,sdata,zarr_name):
  dir_path=os.path.dirname(zarr_name)
  if type!="slide_seq":
    for table in sdata.tables.keys():
      adata = sdata[table]
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
  adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
  adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
  sc.pp.calculate_qc_metrics(
      adata, 
      qc_vars=["mt", "ribo", "hb"], 
      percent_top=(10, 20, 50, 150),
      inplace=True, 
      log1p=True)
  
  if type=="xenium":
    cprobes = (
      adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    cwords = (adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    # print(f"Negative DNA probe count % : {cprobes}")
    # print(f"Negative decoding count % : {cwords}")

  if type=='xenium':
    image_num=4
  else:
    image_num=2
  
  fig, axs = plt.subplots(1, image_num, figsize=(15, 4))
  axs[0].set_title("Total transcripts per cell")
  sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0])

  axs[1].set_title("Unique transcripts per cell")
  sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1])
  if type=='xenium':
    axs[2].set_title("Area of segmented cells")
    sns.histplot(
      adata.obs["cell_area"],
      kde=False,
      ax=axs[2])

    axs[3].set_title("Nucleus ratio")
    sns.histplot(
      adata.obs["nucleus_area"] / adata.obs["cell_area"],
      kde=False,
      ax=axs[3])

  plt.savefig(
      os.path.join(dir_path, "total.png"),
      dpi=300,
      bbox_inches='tight')
  plt.close()
  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts"], 
    groupby="sample", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total UMI by Sample")
  plt.axhline(y=4, color='r', linestyle='-')
  plt.axhline(y=8, color='r', linestyle='-')
  plt.savefig(
    os.path.join(dir_path, "total_umi_by_sample.png"),
    dpi=300, 
    bbox_inches='tight')
  plt.close()



  sc.pl.violin(
    adata=adata, 
    keys=["log1p_n_genes_by_counts"], 
    groupby="sample", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total Genes by Sample")
  plt.savefig(
    os.path.join(dir_path, "total_genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()

  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts_mt"], 
    groupby="sample", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Mitochondrial Genes by Sample")
  plt.savefig(
    os.path.join(dir_path, "genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
  
  
  sc.pl.scatter(adata, "log1p_total_counts_mt", "log1p_n_genes_by_counts", color="pct_counts_mt")
  plt.savefig(
    os.path.join(dir_path, "scatter.png"),
    dpi=300,
    bbox_inches='tight')
  plt.close()
  if type!='slide_seq':
    adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
    adata.obs['region'] = adata.obs['region'].astype('category')
    for table in sdata.tables.keys():
      sdata[table]=adata
  return sdata



sdatas = []
if type=='visium_segment':
  for i in range(len(sample)):
    sdata=spd.read_zarr(args.input_path[i])
    for table in sdata.tables.values():
      table.var_names_make_unique()
      table.obs['cell_id'] = table.obs.index
      table.obs["sample"]=sample[i]
      table.obs["group"]=group[i]
    sdatas.append(sdata)
    del sdata,table
elif type=='visium_HD':
  for i in range(len(sample)):
    sdata=spd.read_zarr(args.input_path[i])
    TABLE_KEY = 'segmentation_counts'
    for shapes in sdata.shapes.keys():
      shape_key=shapes
    for table in sdata.tables.values():
        # table.obs['cell_id'] = table.obs['location_id']
        # table.var_names_make_unique()
        table.obs['cell_id'] = table.obs.index
        table.obs["sample"]=sample[i]
        table.obs["group"]=group[i]
        table.obs['region']=shape_key
        sdata.shapes[shape_key].index=table.obs.index
    del table.uns['spatialdata_attrs']
    sdata.tables={
            TABLE_KEY: TableModel.parse(
                table,
                region=shape_key, # Link table to shapes element
                region_key='region', # Column in adata.obs indicating region name
                instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
            )
        }
    sdatas.append(sdata)
elif type=="visium":
  for i in range(len(sample)):
      sdata=spd.read_zarr(args.input_path[i])
      SHAPES_KEY = sample[i]
      TABLE_KEY = 'segmentation_counts'
      for table in sdata.tables.values():
          table.obs["sample"] = sample[i]
          table.obs["group"]=group[i]
          table.obs['cell_id'] = table.obs.index
          sdata.shapes[sample[i]].index=table.obs['cell_id']
      del table.uns['spatialdata_attrs']
      sdata.tables={
              TABLE_KEY: TableModel.parse(
                  table,
                  region=SHAPES_KEY, # Link table to shapes element
                  region_key='region', # Column in adata.obs indicating region name
                  instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
              )
          }
      sdatas.append(sdata)
elif type=="slide_seq":
  for i in range(len(sample)):
    sdata = sc.read_h5ad(args.input_path[i])
    sdata.obs['cell_id'] = table.obs.index
    sdata.obs['sample'] = sample[i]
    sdata.obs["group"]=group[i]
    sdata.var_names_make_unique()
    sdata.obs_names_make_unique()
    sdatas.append(sdata)
  adata = anndata.concat(sdatas, join='inner', index_unique=None)
  adata.obs_names_make_unique
  adata=QC_plot(type,adata,zarr_name)
  adata.write("./concatenated_sdata")
  exit()
elif type=="xenium":
  for i in range(len(sample)):
    sdata=spd.read_zarr(args.input_path[i])
    new_images = {}
    for img_name in sdata.images.keys():
        new_name = f"{sample[i]}_{img_name}"
        new_images[new_name] = sdata.images[img_name]
    sdata.images = new_images
    new_images = {}
    for shapes_name in sdata.shapes.keys():
        new_name = f"{sample[i]}_{shapes_name}"
        new_images[new_name] = sdata.shapes[shapes_name]
    sdata.shapes = new_images
    new_images = {}
    for labels_name in sdata.labels.keys():
        new_name = f"{sample[i]}_{labels_name}"
        new_images[new_name] = sdata.labels[labels_name]
    sdata.labels = new_images

    new_images = {}
    for points_name in sdata.points.keys():
        new_name = f"{sample[i]}_{points_name}"
        new_images[new_name] = sdata.points[points_name]
    sdata.points = new_images

    new_images = {}
    SHAPES_KEY = f"{sample[i]}_{shapes_name}"
    TABLE_KEY = 'segmentation_counts'
    for table in sdata.tables.values():
        table.var_names_make_unique()
        table.obs["sample"] = sample[i]
        table.obs["group"]=group[i]
        table.obs['cell_id'] = table.obs.index.astype(str)
        table.obs['region'] = SHAPES_KEY
        table.obs['region'] = table.obs['region'].astype('category')
    if SHAPES_KEY in sdata.shapes:
        sdata.shapes[SHAPES_KEY].index = table.obs['cell_id']
    del table.uns['spatialdata_attrs']
    sdata.tables={
            TABLE_KEY: TableModel.parse(
                table,
                region=SHAPES_KEY,
                region_key='region',
                instance_key='cell_id'
            )
        }
    sdatas.append(sdata)

concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
concatenated_sdata=QC_plot(type,concatenated_sdata,args.output_zarr_path)
concatenated_sdata.write(args.output_zarr_path, overwrite=True)
del concatenated_sdata, sdatas
gc.collect()
