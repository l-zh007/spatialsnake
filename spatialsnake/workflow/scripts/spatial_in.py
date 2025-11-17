import os
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
import seaborn as sns

parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--count_file', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=False,
                   help='Path for the output zarr file')
parser.add_argument('--channel', type=str, required=False,
                   help='Path for the output zarr file')                   
########## args of visium HD
parser.add_argument("--bin_size", type=str, required=False, help="bin")

####### args of xenium
parser.add_argument("--cells_boundaries", type=bool, required=False, help="bin")
parser.add_argument("--nucleus_boundaries", type=bool, required=False, help="bin")
parser.add_argument("--nucleus_labels", type=bool, required=False, help="bin")
parser.add_argument("--morphology_mip", type=bool, required=False, help="bin")

##########args of segment visium
parser.add_argument("--scale_factors", type=str, required=False, help="bin")
parser.add_argument("--image", type=str, required=False, help="bin")
parser.add_argument("--geojson", type=str, required=False, help="bin")


######args of slide seq
parser.add_argument("--coor_file", type=str, required=False, help="bin")

#######args of Merfish

args = parser.parse_args()


count_file= args.count_file
real_dir = "/".join(os.path.normpath(args.input_dir).split(os.path.sep)[:2])
type=args.type


def QC_plot(type,sdata,zarr_name):
  dir_path=os.path.dirname(zarr_name)
  if type!="slide_seq":
    for table in sdata.tables.keys():
      adata = sdata[table]
  else:
    adata=sdata
  adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
  sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars="mt",
    percent_top=(10, 20, 50),
    inplace=True, 
    log1p=True)
  
  if type=="xenium":
    cprobes = (
      adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    cwords = (adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)
    print(f"Negative DNA probe count % : {cprobes}")
    print(f"Negative decoding count % : {cwords}")

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
    groupby="region", 
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
  plt.show()
  plt.close()



  sc.pl.violin(
    adata=adata, 
    keys=["log1p_n_genes_by_counts"], 
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Total Genes by Sample")
  plt.savefig(
    os.path.join(dir_path, "total_genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()

  sc.pl.violin(
    adata=adata, 
    keys=["log1p_total_counts_mt"], 
    groupby="region", 
    stripplot=False, 
    inner="box",
    show=False)
  plt.title("Mitochondrial Genes by Sample")
  plt.savefig(
    os.path.join(dir_path, "genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()
  
  
  sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
  plt.savefig(
    os.path.join(dir_path, "scatter.png"),
    dpi=300,
    bbox_inches='tight')
  plt.show()
  plt.close()
  if type!="slide_seq":
    sdata[table]=adata
    return sdata
  else:
    return adata












def create_zarr_bin(path_to_inputs,sample_id,zarr_name,filtered_counts_file,bin_size):
    print(zarr_name)
    sdata = spatialdata_io.visium_hd(path_to_inputs, dataset_id=sample_id, filtered_counts_file=filtered_counts_file, bin_size=bin_size)
    for table in sdata.tables.values():
        table.obs['cell_id'] = table.obs.index
        table.obs["group"] = sample_id
    sdata=QC_plot(type,sdata,zarr_name)
    sdata.write(zarr_name, overwrite=True)

def creat_zarr_visium(path_to_inputs,sample_id,zarr_name,h5_name):
    sdata=spatialdata_io.visium(path_to_inputs,dataset_id=sample_id,counts_file=h5_name)
    sdata=QC_plot(type,sdata,zarr_name)
    SHAPES_KEY = sample_id
    TABLE_KEY = 'table'
    for table in sdata.tables.values():
          table.obs["sample"] = sample_id
          table.obs["group"] = sample_id
          table.obs['cell_id'] = table.obs.index
          sdata.shapes[sample_id].index=table.obs['cell_id']
    del table.uns['spatialdata_attrs']
    sdata.tables={
              TABLE_KEY: TableModel.parse(
                  table,
                  region=SHAPES_KEY, # Link table to shapes element
                  region_key='region', # Column in adata.obs indicating region name
                  instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
              )
          }
    print(sdata["table"].obs)
    sdata.write(zarr_name, overwrite=True)

def create_zarr(count_matrix_path,image_path,scale_factors_path,geojson_path,sample_name,zarr_name):
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
    adata.obs['group'] = sample_name
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
    if args.channel=="single_analysis":
      sdata=QC_plot(type,sdata,zarr_name)
    sdata.write(zarr_name, overwrite=True)

def create_zarr_xenium(path_to_inputs,zarr_name,cells_boundaries,nucleus_boundaries,nucleus_labels,morphology_mip):
    sdata = spatialdata_io.xenium(path_to_inputs,cells_boundaries=cells_boundaries, nucleus_boundaries=nucleus_boundaries,nucleus_labels=nucleus_labels,morphology_mip=morphology_mip,n_jobs=4,cells_as_circles=True)
    for table in sdata.tables.values():
        table.obs['cell_id'] = table.obs.index
        table.obs["group"] = args.sample_id
    sdata=QC_plot(type,sdata,zarr_name)
    sdata.write(zarr_name,overwrite=True)
    
def creat_zarr_slide_seq(zarr_name,count_file,coor_file):
  counts = pd.read_csv(count_file, sep='\t', index_col=0,comment='#')
  coor_df = pd.read_csv(coor_file, index_col=0)
  adata = sc.AnnData(counts.T)
  adata.var_names_make_unique()
  coor_df = coor_df.loc[adata.obs_names, ['xcoord', 'ycoord']]
  adata.obsm["spatial"] = coor_df.to_numpy()
  adata.obs['region']=args.sample_id
  adata.obs['cell_id'] = adata.obs.index
  adata.obs["group"] = args.sample_id
  adata=QC_plot(type,adata,zarr_name)
  adata.write(zarr_name)

if type=='visium_segment':
        count_file=os.path.join(f'data/{args.sample_id}/segmented_outputs',count_file)
        image_file=os.path.join(f'data/{args.sample_id}/segmented_outputs/spatial',args.image)
        scale_factors_file=os.path.join(f'data/{args.sample_id}/segmented_outputs/spatial',args.scale_factors)
        geojson_file=os.path.join(f'data/{args.sample_id}/segmented_outputs',args.geojson)
        print(count_file,args.image,args.scale_factors,args.geojson)
        create_zarr(count_matrix_path=count_file,
                image_path=image_file,
                scale_factors_path=scale_factors_file,
                geojson_path=geojson_file,
                sample_name=args.sample_id,
                zarr_name=args.output_zarr_path)
elif type=='visium_HD':
      is_filtered = True if count_file == "filtered_feature_bc_matrix.h5" else False
      create_zarr_bin(path_to_inputs=real_dir,
                sample_id=args.sample_id,
                zarr_name=args.output_zarr_path,
                filtered_counts_file=is_filtered,
                bin_size=args.bin_size)
elif type=="visium":
      creat_zarr_visium(
                path_to_inputs=real_dir,
                sample_id=args.sample_id,
                zarr_name=args.output_zarr_path,
                h5_name=count_file)
elif type=="slide_seq":
      coor_file=os.path.join(f'data/{args.sample_id}',args.coor_file)
      count_file=os.path.join(f'data/{args.sample_id}',args.count_file)
      creat_zarr_slide_seq(zarr_name=args.output_zarr_path,count_file=count_file,coor_file=coor_file)
elif type=="xenium":
      print(args.cells_boundaries,args.nucleus_boundaries,args.nucleus_labels,args.morphology_mip)
      create_zarr_xenium(path_to_inputs=real_dir,zarr_name=args.output_zarr_path,cells_boundaries=args.cells_boundaries,nucleus_boundaries=args.nucleus_boundaries,nucleus_labels=args.nucleus_labels,morphology_mip=args.morphology_mip)



      









