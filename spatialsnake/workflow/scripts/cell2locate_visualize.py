import os
import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import scanpy as sc
import matplotlib.pyplot as plt
from spatialdata.transformations import Identity, Scale
import argparse
import matplotlib as mpl
from spatialdata.datasets import blobs_annotating_element
from spatialdata.transformations import Affine, set_transformation, get_transformation
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata
import numpy as np
from spatialdata import bounding_box_query
from spatialsnake.workflow.function.transform import zarr_to_h5ad
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
type = args.type

def parse_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    value_str = str(value).strip().lower()
    if value_str in {"true", "1", "yes", "y", "t"}:
        return True
    if value_str in {"false", "0", "no", "n", "f", ""}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

def crop0(x, crs, x1, x2, y1, y2):
    return bounding_box_query(
        x,
        min_coordinate=[x1, y1],
        max_coordinate=[x2, y2],
        axes=("x", "y"),
        target_coordinate_system=crs
    )

def normalize_factor_names(factor_names):
    return [name.replace("+", "plus") for name in factor_names]

def ensure_factor_columns(adata, factor_names):
    normalized_names = normalize_factor_names(factor_names)
    renamed_columns = adata.obs.columns.str.replace("+", "plus", regex=False)
    adata.obs.columns = renamed_columns
    existing_cols = [col for col in normalized_names if col in adata.obs.columns]
    return normalized_names, existing_cols

def plot_abundance_stats(adata, factor_names, save_dir="."):
    """
    Generate stacked bar plots of cell abundance per cluster.
    """
    print("Generating cell abundance stats plot...")
    
    # Check if clusters exist, if not, compute them based on abundance
    if "clusters" not in adata.obs:
            if "q05_cell_abundance_w_sf" in adata.obsm:
                sc.pp.neighbors(adata, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)
                print(f"Computed neighbors using q05_cell_abundance_w_sf with n_neighbors")
                sc.tl.leiden(adata, key_added="clusters")
                sc.tl.umap(adata, min_dist = 0.3, spread = 1)
                sc.pl.umap(adata, color=['clusters'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20)
                plt.savefig(os.path.join(save_dir, "umap.png"), dpi=300, bbox_inches='tight')
                plt.close()
                print(f"Saved UMAP plot to {os.path.join(save_dir, 'umap.png')}")
    valid_factors = [f for f in factor_names if f in adata.obs.columns]
    if not valid_factors:
        print("No valid cell type columns found for plotting.")
        return
    try:
        df = adata.obs[["clusters"] + valid_factors].copy()
        
        # Group by cluster and calculate mean abundance
        cluster_means = df.groupby("clusters")[valid_factors].mean()
        
        # Calculate percentages
        cluster_percents = cluster_means.div(cluster_means.sum(axis=1), axis=0) * 100
        
        # Plotting
        # Setup figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot stacked bar
        cluster_percents.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
        
        ax.set_title("Cell Type Abundance by Cluster")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Percentage Abundance (%)")
        
        # Legend
        ax.legend(title="Cell Types", bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        
        # Save plot
        save_path = os.path.join(save_dir, "cluster_abundance_stacked_bar.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved abundance plot to {save_path}")
        
        # Save stats
        stats_path = os.path.join(save_dir, "cluster_abundance_stats.csv")
        cluster_percents.to_csv(stats_path)
        print(f"Saved abundance stats to {stats_path}")
        
    except Exception as e:
        print(f"Error plotting abundance stats: {e}")
        import traceback
        traceback.print_exc()

def render_spatial_plots(concatenated_sdata, shapes, images, systems, output_dir, image_slice, coord):
    x1, x2, y1, y2 = coord if coord else (None, None, None, None)
    for i in range(len(images)):
        sys = systems[0] if len(systems) < 2 else systems[i]
        axes = plt.subplots(2, 1, figsize=(20, 13))[1].flatten()
        if image_slice and coord:
            concatenated_sdata = crop0(concatenated_sdata, sys, x1, x2, y1, y2)
        concatenated_sdata.pl.render_images(images[i]).pl.show(ax=axes[0], title="image", coordinate_systems=sys)
        concatenated_sdata.pl.render_images(images[i]).pl.render_shapes(shapes[i], color="cellLoca_type").pl.show(
            ax=axes[1], coordinate_systems=sys
        )
        plt.savefig(f"{output_dir}/figure/{images[i]}_Most_rank_celltype.png", dpi=300, bbox_inches="tight")
        plt.close()

image_slice = parse_bool(args.image_slice)
output_dir = os.path.dirname(args.output_zarr_path)

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "figure"), exist_ok=True)

concatenated_sdata = spd.read_zarr(args.input_dir)
table_key = next(iter(concatenated_sdata.tables.keys()))
adata_vis = concatenated_sdata[table_key]
flag = True
adata = zarr_to_h5ad([args.input_dir],f"{output_dir}/test.h5ad",save_image=True)


if "mod" in concatenated_sdata[table_key].uns and "factor_names" in concatenated_sdata[table_key].uns["mod"]:
    factor_names = concatenated_sdata[table_key].uns["mod"]["factor_names"]
    if "q05_cell_abundance_w_sf" in concatenated_sdata[table_key].obsm:
        concatenated_sdata[table_key].obs[factor_names] = concatenated_sdata[table_key].obsm["q05_cell_abundance_w_sf"]
else:
    factor_names = []

normalized_names, existing_cell_cols = ensure_factor_columns(concatenated_sdata[table_key], factor_names)
if len(existing_cell_cols) > 0:
    abundance_matrix = concatenated_sdata[table_key].obs[existing_cell_cols]
    concatenated_sdata[table_key].obs["cellLoca_type"] = abundance_matrix.idxmax(axis=1)
else:
    concatenated_sdata[table_key].obs["cellLoca_type"] = "unknown"

coord = tuple(args.coord) if args.coord else None
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
if shape_count>1 and type!="visium" and args.shape_type:
  for i in range(len(shape_elements)):
    if args.shape_type in shape_elements[i]:
      shapes.append(shape_elements[i])
else:
  shapes=shape_elements

print(shapes)
print(images)
print(systems)
if len(images) != len(shapes):
    print("Check the spatial data to make sure that for every image there is a shape")
    print(images,shapes)
    exit()

else:
    render_spatial_plots(
        concatenated_sdata,
        shapes,
        images,
        systems,
        output_dir,
        image_slice,
        coord
    )

print(adata.obs["sample"])

from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(11, 13), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 16 # number of training restarts
    },
    # the hyperparameters of NMF can be also adjusted:
    model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
    export_args={'path': f'{output_dir}/CoLocatedComb/'}
)





concatenated_sdata.write(output_dir +f"/{args.sample_id}.zarr")
print(concatenated_sdata["Human_Lymph_Node_hires_image"])
print('@@@@@@@@@@@@@@@')
############# each type 可视化
from cell2location.plt import plot_spatial
if flag and len(existing_cell_cols) > 0:
    # 切换当前工作目录到 output_dir/figure 以适应 scanpy 的 save 参数行为
    current_dir = os.getcwd()
    target_fig_dir = os.path.join(output_dir, "figure")
    os.makedirs(target_fig_dir, exist_ok=True)
    os.chdir(target_fig_dir)
    
    try:
        concatenated_sdata[table_key].uns["spatial"]=adata.uns["spatial"]
        library_id = None
        for key in adata.uns["spatial"].keys():
            if "images" in adata.uns["spatial"][key]:
                library_id = key
                break
        if 'scalefactors' in concatenated_sdata[table_key].uns['spatial'][library_id]:
                sf = concatenated_sdata[table_key].uns['spatial'][library_id]['scalefactors']
                # 如果 hires scale 为 1.0 或缺失，尝试重新计算
        if sf.get('tissue_hires_scalef', 1.0) == 1.0:
            print("检测到 tissue_hires_scalef 可能不正确 (1.0)，尝试自动修复...")  
            img = concatenated_sdata[table_key].uns['spatial'][library_id]['images']['hires']
            img_h, img_w = img.shape[:2] if len(img.shape) >= 2 else (0, 0)
            spatial_coords = concatenated_sdata[table_key].obsm['spatial']
            spot_max_x = spatial_coords[:, 0].max()
            spot_max_y = spatial_coords[:, 1].max()
        if img_w > 0 and spot_max_x > img_w:
                scale_x = img_w / spot_max_x
                scale_y = img_h / spot_max_y
                new_scale = min(scale_x, scale_y)-0.025
                sf['tissue_hires_scalef'] = new_scale
                print(f"  -> 启发式估算 scale: {new_scale}")

        concatenated_sdata[table_key].obs[concatenated_sdata[table_key].uns["mod"]["factor_names"]] = concatenated_sdata[table_key].obsm["q05_cell_abundance_w_sf"]
        
        sc.pl.spatial(
            concatenated_sdata[table_key],
            color=existing_cell_cols,
            show=False,
            img_key="hires",
            library_id=library_id,
            size=6
        )
        plt.savefig("each_celltype.png", dpi=300, bbox_inches='tight')
        plt.close()

        sc.pl.spatial(
            concatenated_sdata[table_key],
            color=concatenated_sdata[table_key].uns["mod"]["factor_names"],
            show=False,
            img_key="hires",
            library_id=library_id,
            size=6
        )
        plt.savefig("factor_namescelltype.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Call abundance stats plot
        if "mod" in concatenated_sdata[table_key].uns:
             plot_abundance_stats(concatenated_sdata[table_key], concatenated_sdata[table_key].uns["mod"]["factor_names"], ".")
        
        clust_labels = existing_cell_cols[:6] if len(existing_cell_cols) > 7 else existing_cell_cols
        clust_col = ["".join(str(i)) for i in clust_labels]
        print(list(adata.uns["spatial"].values())[0])
    finally:
        os.chdir(current_dir)