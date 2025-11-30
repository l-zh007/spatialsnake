import os
import sys
import argparse
import spatialdata as spd
import scanpy as sc
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
from spatialsnake.workflow.function.legacy_anndata import transform_to_zarr
import spatialdata_plot as splt







def zarr_to_h5ad(INPUT_list, output_h5ad, save_image=False):
    for file_path in INPUT_list:
        sdata = spd.read_zarr(file_path)
        print(sdata)
    systems=[]
    shape_elements = list(sdata.shapes.keys())
    valid_coord_systems = sorted(sdata.coordinate_systems)
    if len(valid_coord_systems)>1:
      for i in range(len(valid_coord_systems)):
          if valid_coord_systems[i] in shape_elements:
            systems.append(valid_coord_systems[i])
    else:
        systems=valid_coord_systems
    print(systems)
    adata=to_legacy_anndata(sdata, include_images=save_image, coordinate_system=systems[0])
    print(adata)
    sc.pl.spatial(
    adata, 
    color='celltype',
    library_id="Lesional_1_hires_image",
    basis='spatial',
    alpha_img=0.8,
    spot_size=15,
    title="Spatial Plot with Image")
    plt.savefig(
    os.path.join("genes_by_sample.png"),
    dpi=300,
    bbox_inches='tight')
    plt.close()
    adata.write_h5ad(output_h5ad)
    return adata


def h5ad_to_zarr(INPUT_list, output_zarr, base_name):
    for file_path in INPUT_list:
      adata = sc.read_h5ad(file_path)
    sdata = transform_to_zarr(adata,base_name)
    print(sdata)
    sdata.pl.render_images("Lesional_1_hires_image").pl.render_shapes(color="celltype").pl.show()
    plt.savefig(
      os.path.join("genes_by_sample.png"),
      dpi=300,
      bbox_inches='tight')
    plt.close()
    sdata.write(output_zarr)
    return sdata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SpatialData transform')
    parser.add_argument('--INPUT', nargs='+', required=True, 
                       help='INPUT zarr')
    parser.add_argument('--output_dir', type=str, required=True, help='output')
    parser.add_argument('--save_image', type=str ,help='')
    parser.add_argument('--transform_to', type=str, 
                       required=True, help='')
    parser.add_argument('--transform_from', type=str, 
                       required=True, help='')
    
    args = parser.parse_args()
    print(args.INPUT)
    save_image = args.save_image=="True"
    if args.transform_from == args.transform_to:
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(args.INPUT[0]))[0]
    print(base_name)
    if args.transform_to == 'h5ad':
        output_path = os.path.join(args.output_dir, f"{base_name}.h5ad")
        print(output_path)
        zarr_to_h5ad(args.INPUT, output_path, save_image)
    elif args.transform_to == 'zarr':
        output_path = os.path.join(args.output_dir, f"{base_name}.zarr")
        print("to zarr")
        h5ad_to_zarr(args.INPUT, output_path, base_name)



