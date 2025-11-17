import os
import spatialdata as spd
import spatialdata_plot
import spatialdata
import scanpy as sc
import scanpy.external as sce
import json
import gc
import matplotlib.pyplot as plt
import argparse
import anndata
import warnings
warnings.filterwarnings("ignore")
import spatialdata_io
import geosketch as sketch
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon




parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--INPUT_FIlE', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_dir', type=str, required=True,
                   help='Path for the output zarr file')

parser.add_argument('--split_by', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--max_x', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--min_x', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--max_y', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--min_y', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--shape_elements', type=str, required=True,
                   help='Path for the output zarr file')                   
args = parser.parse_args()

print(args.max_y)
print(args.min_x)
print(args.min_y)

# cropped_sdata = sdata_ST8059050.query.bounding_box(
#     min_coordinate=[bb_xmin, bb_ymin],
#     max_coordinate=[bb_xmax, bb_ymax],
#     target_coordinate_system="ST8059050",
# )

def crop0(x,min_x,max_x,min_y,max_y,src):
    return spd.bounding_box_query(
        x,
        min_coordinate=[min_x, min_y],
        max_coordinate=[max_x, max_y],
        axes=("x", "y"),
        target_coordinate_system=src)
    
    
def slice_by_cluster(adata,file_type,barcode,dir_path,concatenated_sdata=None):
  value=adata.obs[barcode].unique()
  for val in value:
    subset = adata.obs[barcode] == val
    subset_adata = adata[subset, :]
    # print(list(subset_adata.obs["region"].unique()))
    # print(subset_adata.uns['spatialdata_attrs'])
    if file_type==".zarr":
      subset_adata.uns['spatialdata_attrs'] = {
    'region': list(subset_adata.obs["region"].unique()),
    'region_key': subset_adata.uns['spatialdata_attrs']["region_key"],
    'instance_key': subset_adata.uns['spatialdata_attrs']["instance_key"]}
      concatenated_sdata[table]=subset_adata
      concatenated_sdata.write(os.path.join(dir_path,f"cluster_{val}.zarr"),overwrite=True)
    else:
      subset_adata.write(os.path.join(dir_path,f"cluster_{val}.h5ad"))
def check_barcode_exit(adata,barcode):
  if barcode in adata.obs.columns:
    return True
  return False

def slice_by_sample(file_type,barcode,concatenated_sdata):
  systems=[]
  shape_elements = list(concatenated_sdata.shapes.keys())
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  if len(valid_coord_systems)>1:
      for i in range(len(valid_coord_systems)):
          if valid_coord_systems[i] in shape_elements:
            systems.append(valid_coord_systems[i])
  else:
      systems=valid_coord_systems
  # print(systems)
  for i in range(len(systems)):
    subset_sample=concatenated_sdata.filter_by_coordinate_system(systems[i])
    # print(subset_sample)
    # print(subset_sample[table].obs)
    subset_sample.write(os.path.join(dir_path,f"cluster_{systems[i]}.zarr"))




dir_path=args.output_dir
file_type=os.path.splitext(args.INPUT_FIlE)[1]
if file_type!=".zarr":
  adata = sc.read_h5ad(args.INPUT_FIlE)
else:
  concatenated_sdata = spd.read_zarr(args.INPUT_FIlE)
  for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]

# print(concatenated_sdata)
# print(adata.obs)
if check_barcode_exit(adata,args.split_by) and args.split_by in ["celltype","cluster","clusters","cell_type","celltypes","cell_types",'clusters','layer','layers']:
  if file_type==".zarr":
    slice_by_cluster(adata,file_type,args.split_by,dir_path,concatenated_sdata=concatenated_sdata)
  else:
    slice_by_cluster(adata,file_type,args.split_by,dir_path)
elif check_barcode_exit(adata,args.split_by) and args.split_by in ["samples","sample","region","group"]:
  if file_type==".zarr":
      slice_by_sample(file_type,args.split_by,concatenated_sdata)
  else:
      slice_by_cluster(adata,file_type,args.split_by,dir_path)
elif args.split_by in ["image","images"]:
  shape_elements = list(concatenated_sdata.shapes.keys())
  valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
  if args.shape_elements and args.shape_elements in shape_elements and args.shape_elements in valid_coord_systems:
      src=args.shape_elements
  else:
      src=valid_coord_systems[0]
  image_id=f"{args.min_x}_{args.max_x}_{args.min_y}_{args.max_y}"
  subset_image=crop0(concatenated_sdata,args.min_x,args.max_x,args.min_y,args.max_y,src)
  subset_image.write(os.path.join(dir_path,f"spatial{image_id}.zarr"),overwrite=True)
  axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
  subset_image.pl.render_images().pl.show(ax=axes[0], title="image",coordinate_systems=src)
  subset_image.pl.render_images().pl.render_shapes(color="clusters").pl.show(ax=axes[1],coordinate_systems=src)
  # subset_image.pl.render_shapes(color="clusters").pl.show()
  plt.savefig(
        os.path.join(dir_path, f"{image_id}_shape.png"),
        dpi=300,
        bbox_inches='tight')
  plt.close()



# sample_cnt=len(adata.obs[barcode].unique())
#   image_elements = list(concatenated_sdata.images.keys())
#   shape_elements = list(concatenated_sdata.shapes.keys())
#   valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
#   image_counts=len(image_elements)//sample_cnt
#   shape_count=len(shape_elements)//sample_cnt
#   shapes=[]
#   images=[]
#   systems=[]
#   if len(valid_coord_systems)>1:
#       for i in range(len(image_elements)):
#         if args.image_type in image_elements[i]:
#           images.append(image_elements[i])
#   else:
#       images=image_elements
# 
#   if len(valid_coord_systems)>1:
#       for i in range(len(valid_coord_systems)):
#           if args.image_type in valid_coord_systems[i]:
#             systems.append(valid_coord_systems[i])
#   else:
#       systems=valid_coord_systems
#   print(shape_count,shape_elements)
#   if shape_count>1 and type!="visium":
#       for i in range(len(shape_elements)):
#         if args.shape_type in shape_elements[i]:
#           shapes.append(shape_elements[i])
#   else:
#       shapes=shape_elements
#   for i in range(len(image_elements)):
#     extent =  spd.get_extent(concatenated_sdata,elements=[shape_elements[i]],coordinate_system='downscale_to_hires')
#     extents.append(extent)
#   for i in range(len(extents)):
#     concatenated_sdata=crop0(concatenated_sdata,crs="Non_Lesional_1_downscaled_hires",bbox=extents[0])
#     concatenated_sdata.write(os.path.join(dir_path,f"sample_{val}.zarr"))





# # axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
# # concatenated_sdata.pl.render_images('Non_Lesional_1_hires_image').pl.show(ax=axes[0], title="image",coordinate_systems='Non_Lesional_1_downscaled_hires')
# # concatenated_sdata.pl.render_images('Non_Lesional_1_hires_image').pl.render_shapes("Non_Lesional_1",color="clusters").pl.show(ax=axes[1],coordinate_systems='Non_Lesional_1_downscaled_hires', title="sd")
# # plt.savefig(os.path.join("tmp_png", f"{val}_Clusters.png"),dpi=300,bbox_inches='tight')
# # plt.close()
#   sample_cnt=args.sample_cnt
#   image_elements = list(concatenated_sdata.images.keys())
#   shape_elements = list(concatenated_sdata.shapes.keys())
#   valid_coord_systems = sorted(concatenated_sdata.coordinate_systems)
#   image_counts=len(image_elements)//sample_cnt
#   shape_count=len(shape_elements)//sample_cnt
#   shapes=[]
#   images=[]
#   systems=[]
#   if len(valid_coord_systems)>1:
#       for i in range(len(image_elements)):
#         if args.image_type in image_elements[i]:
#           images.append(image_elements[i])
#   else:
#       images=image_elements
# 
#   if len(valid_coord_systems)>1:
#       for i in range(len(valid_coord_systems)):
#           if args.image_type in valid_coord_systems[i]:
#             systems.append(valid_coord_systems[i])
#   else:
#       systems=valid_coord_systems
#   print(shape_count,shape_elements)
#   if shape_count>1 and type!="visium":
#       for i in range(len(shape_elements)):
#         if args.shape_type in shape_elements[i]:
#           shapes.append(shape_elements[i])
#   else:
#       shapes=shape_elements









