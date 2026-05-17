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
import re

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon
import argparse
import itertools
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
algorithm = 'edgeR' if args.sample_cnt < 3 else args.algorithm

focus = args.cell_focus
if focus is not None and str(focus).strip().lower() in {"none", "null", ""}:
    focus = None
focus_terms = None
if focus is not None:
    focus_terms = [term.strip() for term in str(focus).split(",") if term.strip() != ""]
    if len(focus_terms) == 0:
        focus_terms = None
print(args.sample_id,args.group)
concatenated_sdata = spd.read_zarr(args.input_dir)
for table in concatenated_sdata.tables.keys():
    table=table
    adata = concatenated_sdata[table]


aggregated = sc.get.aggregate(adata, by=['celltype',"region"], func=["sum"])
print(aggregated.layers["sum"])

print(focus_terms if focus_terms is not None else focus)
metadata_df = pd.DataFrame(
    data={"sample":args.sample_id,"condition": args.group},
    index=args.sample_id
)




if focus_terms:
    celltype_series = aggregated.obs["celltype"].astype(str)
    pattern = "|".join([re.escape(term) for term in focus_terms])
    mask = celltype_series.str.contains(pattern, regex=True)
    aggregated_cluster_of_interest = aggregated[mask]
else:
    aggregated_cluster_of_interest = aggregated
print(metadata_df.index)
print(aggregated.obs.index)

if "region" in aggregated_cluster_of_interest.obs.columns:
    region_series = aggregated_cluster_of_interest.obs["region"].astype(str)
    missing_samples = [s for s in args.sample_id if s not in set(region_series)]
    print(set(args.sample_id), set(region_series))
    if missing_samples:
        print(f"'{args.cell_focus}' not exit in the samples:")
        for sample in missing_samples:
            print(f"  - {sample}")
        print(f"\nplease select a common celltype to compare analyze")
        sys.exit(1)

if "region" in aggregated_cluster_of_interest.obs.columns:
    region_series = aggregated_cluster_of_interest.obs["region"].astype(str)
    counts_df = pd.DataFrame(
        aggregated_cluster_of_interest.layers["sum"],
        index=region_series,
        columns=aggregated_cluster_of_interest.var_names
    )
    counts_df = counts_df.groupby(level=0).sum()
    counts_df = counts_df.reindex(metadata_df.index)
else:
    counts_df = pd.DataFrame(
        aggregated_cluster_of_interest.layers["sum"],
        index=metadata_df.index,
        columns=aggregated_cluster_of_interest.var_names
    )
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
conditions = list(pd.unique(metadata_df["condition"]))
if len(conditions) < 2:
    print("Need at least two conditions for differential expression.")
    sys.exit(1)

contrast_pairs = list(itertools.combinations(conditions, 2))
diff_dir = os.path.join(os.path.dirname(args.output_zarr_path), "diff")
os.makedirs(diff_dir, exist_ok=True)

all_results = []
for group1, group2 in contrast_pairs:
    ds = DeseqStats(dds, contrast=["condition", group1, group2])
    ds.summary()
    result_df = ds.results_df.copy()
    result_df["contrast"] = f"{group1}_vs_{group2}"
    result_df["group1"] = group1
    result_df["group2"] = group2
    result_df["comparison_group"] = group1
    result_df["reference_group"] = group2
    result_df["higher_in_group"] = np.where(
        result_df["log2FoldChange"] > 0,
        group1,
        np.where(result_df["log2FoldChange"] < 0, group2, "None")
    )
    result_df["lower_in_group"] = np.where(
        result_df["log2FoldChange"] > 0,
        group2,
        np.where(result_df["log2FoldChange"] < 0, group1, "None")
    )
    result_df["regulation_label"] = np.where(
        result_df["log2FoldChange"] > 0,
        f"Higher_in_{group1}",
        np.where(result_df["log2FoldChange"] < 0, f"Higher_in_{group2}", "None")
    )
    result_df["gene"] = result_df.index
    result_df.to_csv(os.path.join(diff_dir, f"{group1}_vs_{group2}.csv"))
    all_results.append(result_df)

combined_results = pd.concat(all_results, axis=0)
combined_results.to_csv(args.output_zarr_path)
