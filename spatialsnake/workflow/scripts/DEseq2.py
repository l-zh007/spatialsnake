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
import sys
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
parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
parser.add_argument('--input_dir', type=str, required=True, 
                   help='Path to the raw data directory')
parser.add_argument('--output_zarr_path', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--sample_id', type=str, required=True,nargs='+',
                   help='Path for the output zarr file')
parser.add_argument('--type', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--algorithm', type=str, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--cell_focus', type=str,required=False,
                   help='Path for the output zarr file')
parser.add_argument('--sample_cnt', type=int, required=True,
                   help='Path for the output zarr file')
parser.add_argument('--group', type=str, required=True,nargs='+',
                   help='Path for the output zarr file')                   
              
args = parser.parse_args()
print(args.sample_cnt)
algorithm = 'edgeR' if args.sample_cnt<3 else args.algorithm

focus=args.cell_focus
print(args.sample_id,args.group)
concatenated_sdata = spd.read_zarr(args.input_dir)
for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]


aggregated = sc.get.aggregate(adata, by=['celltype',"region"], func=["sum"])
print(aggregated.layers["sum"])

sample=[f"{args.cell_focus}_{i}" for i in args.sample_id]
print(sample)
metadata_df = pd.DataFrame(
    data={"sample":args.sample_id,"condition": args.group},
    index=args.sample_id
)




if focus:
    aggregated_cluster_of_interest = aggregated[aggregated.obs["celltype"] == args.cell_focus]
else:
    aggregated_cluster_of_interest = aggregated
print(metadata_df.index)
print(aggregated.obs.index)

all_samples = set(sample)
available_samples = set(aggregated_cluster_of_interest.obs.index)
missing_samples = all_samples - available_samples
print(all_samples,available_samples)
if missing_samples:
    print(f"'{args.cell_focus}' not exit in the samples:")
    for sample in missing_samples:
        print(f"  - {sample}")
    print(f"\nplease select a common celltype to compare analyze")
    sys.exit(1)

counts_df = pd.DataFrame(aggregated_cluster_of_interest.layers["sum"],index=metadata_df.index,columns=aggregated_cluster_of_interest.var_names)
counts_df = counts_df.astype(int)
print(counts_df)
print(metadata_df)


if algorithm=="edgeR":
  metadata_df.to_csv(os.path.join(os.path.dirname(args.output_zarr_path),"group.csv"))
  counts_df = counts_df.T
  counts_df.to_csv(
    args.output_zarr_path,
    index=True,
    header=True,
    float_format=None)
  exit()

dds = DeseqDataSet(
    counts = counts_df,
    metadata = metadata_df,
    design = "~condition",
    refit_cooks=True
)

dds.deseq2()
ds = DeseqStats(dds, contrast=["condition", "Cancer", "Normal"])
ds.summary()

print(ds.results_df.sort_values(by='log2FoldChange', ascending=False).head(10))

positive_fold_change_genes = ds.results_df[ds.results_df['log2FoldChange'] > 0]
log2fc_threshold = 1.0
padj_threshold = 0.5

ds.results_df['significant'] = (
    (abs(ds.results_df['log2FoldChange']) > log2fc_threshold) & 
    (ds.results_df['padj'] < padj_threshold)
)

ds.results_df.to_csv(args.output_zarr_path)











# top_10_positive_genes_table = positive_fold_change_genes.sort_values(by='padj', ascending=True).head(10)
# 
# print("\n--- Top 10 Positively Differentially Expressed Genes (Cancer vs. Normal) ---")
# print(top_10_positive_genes_table)
# 
# 
# 
# 
# image_elements = list(concatenated_sdata.images.keys())
# shape_elements = list(concatenated_sdata.shapes.keys())
# gene_name = "COL1A1"
# extents = []
# 
# for i in range(len(image_elements)):
#     extent =  spd.get_extent(concatenated_sdata,elements=[shape_elements[i]],coordinate_system='downscale_to_hires')
#     extents.append(extent)
# for i in range(len(image_elements)):
#   print("Plotting: "+ image_elements[i])
#   title=image_elements[i].replace("_gene_image","")
#   crop0(concatenated_sdata,crs="downscale_to_hires",bbox=extents[i]).pl.render_images(image_elements[i]).pl.render_shapes(shape_elements[i],color=gene_name).pl.show(coordinate_systems="downscale_to_hires",title=title)
#   plt.savefig(
#         os.path.join(output_dir, f"{title}_{gene_name}_express.png"),
#         dpi=300,
#         bbox_inches='tight')
#   plt.close()





# # 获取 AnnData 对象的 obs 数据框
# obs_df = concatenated_sdata["segmentation_counts"].obs
# 
# # 步骤1：添加新类别 "T cell" 到 grouped_clusters 的类别列表中
# obs_df["grouped_clusters"] = obs_df["grouped_clusters"].cat.add_categories("T cell")
# 
# # 步骤2：用 "T cell" 填充缺失值
# obs_df["grouped_clusters"] = obs_df["grouped_clusters"].fillna("T cell")
# 
# 
# concatenated_sdata["segmentation_counts"].obs = obs_df
# 
# print(concatenated_sdata["segmentation_counts"].obs["clusters"][1])





