import os
import sys
import argparse
import spatialdata as spd
import scanpy as sc
import geopandas as gpd
from shapely.geometry import Point
import numpy as np

def zarr_to_h5ad(zarr_path, output_h5ad, save_image=False):
    sdata = spd.read_zarr(zarr_path)
    if not sdata.tables:
        raise ValueError("No tables found in SpatialData object")
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name].copy()
    if hasattr(sdata, 'shapes') and len(sdata.shapes) > 0:
        shape_key = list(sdata.shapes.keys())[0]
        shapes = sdata.shapes[shape_key]
        coords = np.array([[geom.centroid.x, geom.centroid.y] 
                          for geom in shapes.geometry])
        adata.obsm['spatial'] = coords
        adata.uns['spatial'] = {'shapes_key': shape_key}
    adata.write_h5ad(output_h5ad)
    return adata


def h5ad_to_zarr(h5ad_path, output_zarr, save_image=False):
    adata = sc.read_h5ad(h5ad_path)
    sdata = spd.SpatialData()
    sdata.tables['table'] = adata
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
        geometry = [Point(coord[0], coord[1]) for coord in coords]
        gdf = gpd.GeoDataFrame(
            {'cell_id': adata.obs_names},
            geometry=geometry,
            index=adata.obs_names)
        sdata.shapes['cells'] = gdf
    sdata.write(output_zarr)
    return sdata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SpatialData 单文件格式转换')
    parser.add_argument('--INPUT',nargs='+',type=str, required=True, help='输入文件路径')
    parser.add_argument('--output_dir', type=str, required=True, help='输出目录')
    parser.add_argument('--save_image', type=str ,help='保留图像数据')
    parser.add_argument('--transform_to', type=str, 
                       required=True, help='目标格式')
    parser.add_argument('--transform_from', type=str, 
                       required=True, help='源格式')
    
    args = parser.parse_args()
    print(args.INPUT[0])
    save_image = args.save_image=="True"
    if args.transform_from == args.transform_to:
        print("错误: 源格式和目标格式不能相同")
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(args.INPUT[0]))[0]
    print(base_name)
    if args.transform_to == 'h5ad':
        output_path = os.path.join(args.output_dir, f"{base_name}.h5ad")
        zarr_to_h5ad(args.INPUT[0], output_path, save_image)
    else:
        output_path = os.path.join(args.output_dir, f"{base_name}.zarr")
        h5ad_to_zarr(args.INPUT[0], output_path, save_image)































