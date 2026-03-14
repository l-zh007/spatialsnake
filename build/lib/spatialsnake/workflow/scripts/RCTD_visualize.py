import argparse
import os
import pandas as pd
import spatialdata as spd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import warnings
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv

warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser(description="Visualize RCTD results on SpatialData zarr")
    parser.add_argument("--zarr_input", type=str, required=True, help="Path to input SpatialData zarr")
    parser.add_argument("--rctd_results", type=str, required=True, help="Path to RCTD results CSV (weights)")
    parser.add_argument("--output_zarr", type=str, required=True, help="Path to output annotated SpatialData zarr")
    parser.add_argument("--sample_id", type=str, required=True, help="Sample ID")
    parser.add_argument("--output_plot", type=str, required=True, help="Path to save the visualization plot")
    
    args = parser.parse_args()

    # 1. Load SpatialData
    print(f"Loading SpatialData from: {args.zarr_input}")
    sdata = spd.read_zarr(args.zarr_input)
    
    # Assuming standard table structure
    # Try to find the table (AnnData)
    table_keys = list(sdata.tables.keys())
    if not table_keys:
        raise ValueError("No tables found in SpatialData zarr")
    
    table_name = table_keys[0] # Default to the first table
    adata = sdata[table_name]
    print(f"Using table: {table_name}")

    # 2. Load RCTD Results
    print(f"Loading RCTD results from: {args.rctd_results}")
    rctd_df = pd.read_csv(args.rctd_results, index_col=0)
    rctd_df.index = rctd_df.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    
    # 3. Merge Results into AnnData
    # Check for common indices
    common_indices = adata.obs_names.intersection(rctd_df.index)
    if len(common_indices) == 0:
        print("Warning: No common indices found between SpatialData and RCTD results!")
        print("AnnData indices example:", adata.obs_names[:5])
        print("RCTD indices example:", rctd_df.index[:5])
        # Try cleaning indices (e.g. replacing dots with dashes if needed, but usually Seurat <-> SpatialData should match if processed correctly)
        # Attempt to match anyway
    else:
        print(f"Found {len(common_indices)} common spots/cells.")

    if rctd_df.shape[1] == 0:
        raise ValueError("RCTD results contain no cell type columns.")

    rctd_df_reindexed = rctd_df.reindex(adata.obs_names)
    
    # Identify cell type columns (assuming all columns in weights csv are cell types)
    cell_types = rctd_df.columns.tolist()
    
    # Add to obsm (easier to manage multiple cell types)
    rctd_weights = rctd_df_reindexed.fillna(0)
    adata.obsm["RCTD_weights"] = rctd_weights
    
    # Also add top prediction to obs for simple plotting
    # Get the column name with max value for each row
    adata.obs["RCTD_max_celltype"] = rctd_weights.idxmax(axis=1).astype(str)
    adata.obs["RCTD_max_weight"] = rctd_weights.max(axis=1)

    # Add individual cell types to obs for individual plotting if needed (optional, can clutter obs)
    # for ct in cell_types:
    #     adata.obs[f"RCTD_{ct}"] = rctd_df_reindexed[ct]

    # Update sdata table
    # In-memory update, sdata holds reference to adata
    
    # 4. Save updated Zarr
    # Since we are modifying the table, we might need to write it back.
    # If output_zarr is same as input, we overwrite. 
    # SpatialData write support is evolving. 
    # Safest way often is to save the modified AnnData back to the Zarr group or write a new Zarr.
    
    # For this pipeline, usually we write to a new location or overwrite the 'annotion' slot
    print(f"Saving annotated SpatialData to: {args.output_zarr}")
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output_zarr), exist_ok=True)

    export_dir = os.path.dirname(args.output_zarr) or "."
    data_type = "xenium" if "xenium" in f"{args.zarr_input}{args.output_zarr}".lower() else "visium"
    export_cluster_csv(
        sdata,
        data_type,
        export_dir,
        cell_id_col="cell_id",
        info_col="RCTD_max_celltype",
        sample_col="sample",
        sample_id=args.sample_id
    )
    
    # If paths are different, we might want to copy the whole zarr first or write a new one.
    # However, sdata.write() writes the whole object.
    # Note: sdata.write(overwrite=True) might be needed.
    
    # Workaround: If we want to save just the table or the whole sdata structure
    # For simplicity in this pipeline context, we assume we write the whole sdata structure to the new location.
    sdata.write(args.output_zarr, overwrite=True)

    # 5. Visualization
    print(f"Generating plots to: {args.output_plot}")
    
    # 提取可视化所需的元素
    image_elements = list(sdata.images.keys())
    shape_elements = list(sdata.shapes.keys())
    valid_coord_systems = sorted(sdata.coordinate_systems)
    
    # 简单的逻辑：选择第一个可用的图像和形状进行绘图
    # 在真实场景中可能需要更复杂的逻辑来匹配最佳的图像/形状对
    
    shapes_to_plot = []
    images_to_plot = []
    systems_to_plot = []
    
    # 尝试找到 hires 图像，如果没有则使用第一个
    hires_images = [img for img in image_elements if "hires" in img]
    if hires_images:
        images_to_plot = hires_images
    elif image_elements:
        images_to_plot = [image_elements[0]]
        
    # 找到关联的 coordinate system
    if images_to_plot:
        # 假设每个图像对应一个坐标系，这里简单处理
        # 实际上应该检查 sdata.images[img].attrs['transform'] 等信息
        # 但在 spatialdata 库中，通常可以直接通过 coordinate_systems 参数指定
        pass
        
    # 确定要绘制的形状 (通常是 cell_boundaries 或 nucleus_boundaries)
    # 如果没有特定的形状，就用第一个
    if shape_elements:
         shapes_to_plot = [shape_elements[0]]

    # 确定坐标系
    if valid_coord_systems:
        systems_to_plot = [valid_coord_systems[0]]
    
    # 绘图逻辑
    # 如果有图像和形状，叠加绘制
    # 如果只有形状（如 Slide-seq），只绘制形状
    
    if not images_to_plot and not shapes_to_plot:
        print("Warning: No images or shapes found to plot.")
        return

    # 准备绘图
    # 这里我们只画第一组匹配的图像/形状
    
    img_key = images_to_plot[0] if images_to_plot else None
    shape_key = shapes_to_plot[0] if shapes_to_plot else None
    sys_key = systems_to_plot[0] if systems_to_plot else None

    fig, ax = plt.subplots(figsize=(10, 10))
    plotted = False

    if img_key:
        try:
            if sys_key:
                sdata.pl.render_images(img_key).pl.show(ax=ax, coordinate_systems=sys_key, title="")
            else:
                sdata.pl.render_images(img_key).pl.show(ax=ax, title="")
            plotted = True
        except Exception as e:
            print(f"Warning: Failed to render image {img_key}: {e}")

    if shape_key:
        try:
            if sys_key:
                sdata.pl.render_shapes(shape_key, color="RCTD_max_celltype").pl.show(ax=ax, coordinate_systems=sys_key, title=f"{args.sample_id} - RCTD Top Cell Type")
            else:
                sdata.pl.render_shapes(shape_key, color="RCTD_max_celltype").pl.show(ax=ax, title=f"{args.sample_id} - RCTD Top Cell Type")
            plotted = True
        except Exception as e:
            print(f"Warning: Failed to render shapes {shape_key}: {e}")

    if not plotted:
        coords = None
        if "spatial" in adata.obsm:
            coords = adata.obsm["spatial"]
        elif all(c in adata.obs.columns for c in ["x", "y"]):
            coords = adata.obs[["x", "y"]].to_numpy()
        if coords is None:
            print("Warning: No valid coordinate system or coordinates found for plotting.")
            return
        labels = pd.Categorical(adata.obs["RCTD_max_celltype"])
        ax.scatter(coords[:, 0], coords[:, 1], c=labels.codes, s=5, cmap="tab20")
        ax.set_title(f"{args.sample_id} - RCTD Top Cell Type")
        ax.set_aspect("equal")

    plt.savefig(args.output_plot, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Plot saved to {args.output_plot}")

if __name__ == "__main__":
    main()
