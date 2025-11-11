import os
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import gc
import matplotlib.pyplot as plt
import argparse
import anndata
import warnings
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--INPUT_FIlE', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
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
args = parser.parse_args()

print(args.INPUT_FIlE)
print(args.output_zarr_path)
print(args.split_by)

# dir_path=os.path.dirname(args.output_zarr_path)
# file_type=os.path.splitext(args.input_path)
# if file_type==".zarr":
#   adata = sc.read_h5ad(args.input_path)
# else:
#   concatenated_sdata = spd.read_zarr(args.input_path)
#   for table in concatenated_sdata.tables.keys():
#     table=table
#     adata = concatenated_sdata[table]
# 
# def crop0(x,crs,bbox):
#     return spd.bounding_box_query(
#         x,
#         min_coordinate=[bbox['x'][0], bbox['y'][0]],
#         max_coordinate=[bbox['x'][1], bbox['y'][1]],
#         axes=("x", "y"),
#         target_coordinate_system=crs,
#     )
# 
# def slice_by_cluster(file_type,barcode):
#   value=adata.obs[barcode].unique()
#   for val in value:
#     subset = adata.obs[barcode] == val
#     subset_adata = adata[subset, :]
#     if file_type==".zarr":
#       concatenated_sdata[table]=subset_adata
#       concatenated_sdata.write(os.path.join(dir_path,f"cluster_{val}.zarr"))
#     else:
#       subset_adata.write(os.path.join(dir_path,f"cluster_{val}.h5ad"))
# 
# 
# def slice_by_sample(file_type):
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
#   for i in range(len(image_elements)):
#     extent =  spd.get_extent(concatenated_sdata,elements=[shape_elements[i]],coordinate_system='downscale_to_hires')
#     extents.append(extent)
#   for i in range(len(extents)):
#     concatenated_sdata=crop0(concatenated_sdata,crs="Non_Lesional_1_downscaled_hires",bbox=extents[0])
#     concatenated_sdata.write(os.path.join(dir_path,f"sample_{val}.zarr"))
# 
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









