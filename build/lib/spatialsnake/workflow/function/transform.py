import os
import sys
import subprocess
import argparse
import spatialdata as spd
import scanpy as sc
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt
from spatialsnake.workflow.function.legacy_anndata import transform_to_zarr
from spatialsnake.workflow.function.legacy_anndata import to_legacy_anndata
import json
def zarr_to_h5ad(INPUT_list, output_h5ad, save_image=True):
    all_adatas = []

    for file_path in INPUT_list:
        print(f"Reading Zarr: {file_path}")
        sdata = spd.read_zarr(file_path)
        print(sdata)
        systems=[]
        shape_elements = list(sdata.shapes.keys())
        for i in range(len(shape_elements)):         ### HD shape str cleaning 
            if "_square" in shape_elements[i]:
                shape_elements[i]=shape_elements[i].split('_square')[0]
        valid_coord_systems = sorted(sdata.coordinate_systems)
        
        if len(valid_coord_systems)>1:
          for i in range(len(valid_coord_systems)):
            if valid_coord_systems[i] in shape_elements:
                systems.append(valid_coord_systems[i])
        else:
            systems=valid_coord_systems
        
        print(f"Found systems (samples): {systems}")
        print(f"Found valid_coord_systems: {valid_coord_systems}")
        
        for system in systems:
            print(f"Converting system: {system}")
            try:
                adata = to_legacy_anndata(sdata, include_images=save_image, coordinate_system=system)
                if 'spatial' in adata.uns:
                    spatial_dict = adata.uns['spatial']
                    lib_id = system
                    if lib_id not in spatial_dict:
                        spatial_dict[lib_id] = {}
                    if 'metadata' not in spatial_dict[lib_id]:
                        for key in list(spatial_dict.keys()):
                            if key == lib_id:
                                continue
                            if 'metadata' in spatial_dict.get(key, {}):
                                spatial_dict[lib_id]['metadata'] = spatial_dict[key]['metadata']
                                break
                    images_dict = spatial_dict[lib_id].get('images', {})
                    scalefactors_dict = spatial_dict[lib_id].get('scalefactors', {})
                    keys_to_remove = []
                    for key in list(spatial_dict.keys()):
                        if key == lib_id:
                            continue
                        entry = spatial_dict.get(key, {})
                        if 'images' in entry:
                            if 'hires' in entry['images']:
                                images_dict['hires'] = entry['images']['hires']
                            if 'lowres' in entry['images']:
                                images_dict['lowres'] = entry['images']['lowres']
                            keys_to_remove.append(key)
                        if 'scalefactors' in entry:
                            scalefactors_dict = entry['scalefactors']
                            keys_to_remove.append(key)
                    for key in keys_to_remove:
                        if key in spatial_dict and key != lib_id:
                            del spatial_dict[key]
                    for k, v in list(scalefactors_dict.items()):
                        if isinstance(v, (np.floating, np.integer)):
                            scalefactors_dict[k] = v.item()
                    spatial_dict[lib_id]['images'] = images_dict
                    spatial_dict[lib_id]['scalefactors'] = scalefactors_dict
                    adata.uns['spatial'] = spatial_dict
                    adata.obs['library_id'] = lib_id
                    adata.obs['library_id'] = adata.obs['library_id'].astype('category')
                    print(f"  Assigned library_id: {lib_id}")
                print(adata.obsm['spatial'])
                all_adatas.append(adata)
            except Exception as e:
                print(f"  Error converting system {system}: {e}")

    if not all_adatas:
        print("No AnnData objects created.")
        return None

    if len(all_adatas) > 1:
        print(f"Concatenating {len(all_adatas)} samples...")
        final_adata = sc.concat(all_adatas, join='outer', label='batch', keys=None, index_unique=None)
        
        # Merge uns['spatial']
        final_adata.uns['spatial'] = {}
        for ad in all_adatas:
            if 'spatial' in ad.uns:
                final_adata.uns['spatial'].update(ad.uns['spatial'])
    else:
        final_adata = all_adatas[0]
        
    print(final_adata)
    if 'spatial' in final_adata.uns:
        print(final_adata.uns['spatial'].keys())
        lib_id = final_adata.obs["library_id"].iloc[0] if "library_id" in final_adata.obs else None
        if lib_id is not None and lib_id in final_adata.uns['spatial']:
            print(final_adata.uns['spatial'][lib_id].keys())
    if "library_id" in final_adata.obs:
        print(final_adata.obs["library_id"])
    final_adata.write_h5ad(output_h5ad)
    return final_adata


def h5ad_to_seurat(INPUT_list, output_dir, data_type='st'):
    script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "h5ad_to_seurat.R")
    
    for input_file in INPUT_list:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_rds = os.path.join(output_dir, f"{base_name}.rds")
        
        cmd = ["Rscript", script_path, input_file, output_rds, data_type]
        print(f"Running command: {' '.join(cmd)}")
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            print(f"Error converting {input_file} to Seurat: {e}")
            sys.exit(1)


def h5ad_to_zarr(INPUT_list, output_zarr, base_name):
    for file_path in INPUT_list:
      adata = sc.read_h5ad(file_path)
    sdata = transform_to_zarr(adata,base_name)
    print(sdata)
    # sdata.pl.render_images("Lesional_1_hires_image").pl.render_shapes(color="celltype").pl.show()
    # plt.savefig(
    #   os.path.join("genes_by_sample.png"),
    #   dpi=300,
    #   bbox_inches='tight')
    # plt.close()
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
    parser.add_argument('--type', type=str, default='st',
                       help='Data type for Seurat conversion: sc (single-cell) or st (spatial transcriptomics). Default is st.')
    
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
    elif args.transform_to == 'seurat':
        print(f"Converting to Seurat in directory: {args.output_dir}")
        if args.transform_from == 'h5ad':
            h5ad_to_seurat(args.INPUT, args.output_dir, args.type)
        elif args.transform_from == 'zarr':
            h5ad_output_path = os.path.join(args.output_dir, f"{base_name}.h5ad")
            print(f"Step 1: Converting Zarr to intermediate integrated H5AD: {h5ad_output_path}")
            zarr_to_h5ad(args.INPUT, h5ad_output_path, save_image)
            
            print(f"Step 2: Converting intermediate H5AD to Seurat RDS...")
            if os.path.exists(h5ad_output_path):
                h5ad_to_seurat([h5ad_output_path], args.output_dir, args.type)
            else:
                print(f"Error: Intermediate H5AD file {h5ad_output_path} was not created.")
                sys.exit(1)
