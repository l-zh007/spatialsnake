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
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt
from spatialdata.transformations import Identity, Scale
import argparse
# parser = argparse.ArgumentParser(description='Process spatial data and convert to zarr format')
# parser.add_argument('--INPUT', nargs='+', required=True, 
#                    help='Path to the raw data directory')
# parser.add_argument('--output_dir', type=str, required=True,
#                    help='Path for the output zarr file')
# parser.add_argument('--merge_by', type=str, required=True,
#                    help='Path for the output zarr file')
# parser.add_argument('--reordering', type=str, required=True,
#                    help='Path for the output zarr file')
# parser.add_argument('--re_sample', type=str, required=True,
#                    help='Path for the output zarr file')
# args = parser.parse_args()

def add_sample_barcode(adata,file_path):
  if 'sample' not in adata.obs.columns:
      file_name = os.path.basename(file_path)
      sample_barcode = os.path.splitext(file_name)[0]
      adata.obs['sample'] = sample_barcode
  return adata


def reorder_by_cluster(table, cluster_key, site_counter):
    if cluster_key not in table.obs.columns:
        return False, site_counter, None
    raw_order = table.obs[cluster_key].unique()
    mapping = {val: idx + site_counter for idx, val in enumerate(raw_order)}
    table.obs['clusters'] = table.obs[cluster_key].map(mapping)
    return True, site_counter + len(raw_order), table


def merge_by_samples(INPUT_list, re_sample):
    sdatas = []
    for file_path in INPUT_list:
        sdata = spd.read_zarr(file_path)
        TABLE_KEY = list(sdata.tables.keys())[0]
        table = sdata.tables[TABLE_KEY]
        if re_sample:
            table = add_sample_barcode(table, file_path)
        table.obs['cell_id'] = table.obs.index
        table.uns['spatialdata_attrs']['instance_key'] = 'cell_id'
        sdata.tables[TABLE_KEY] = table
        sdatas.append(sdata)
    return sdatas


def merge_by_clusters(INPUT_list, reordering, cluster_key):
    sdatas = []
    site_counter = 0
    for file_path in INPUT_list:
        sdata = spd.read_zarr(file_path)
        TABLE_KEY = list(sdata.tables.keys())[0]
        table = sdata.tables[TABLE_KEY]
        print(table)
        if reordering and cluster_key in table.obs.columns:
            success, site_counter, updated_table = reorder_by_cluster(
                table, cluster_key, site_counter
            )
            if not success:
                print(f"'{cluster_key}' not in the table")
                sys.exit(1)
            table = updated_table
        table.obs['cell_id'] = table.obs.index
        table.uns['spatialdata_attrs']['instance_key'] = 'cell_id'
        sdata.tables[TABLE_KEY] = table
        sdatas.append(sdata)
    return sdatas


def parse_csv_inputs(annotation_csv):
    raw = str(annotation_csv).strip()
    if raw == "":
        return []
    parts = [p.strip() for p in raw.split(",") if p.strip() != ""]
    paths = []
    for part in parts:
        if os.path.isdir(part):
            csvs = [os.path.join(part, name) for name in os.listdir(part) if name.lower().endswith(".csv")]
            csvs.sort()
            paths.extend(csvs)
        else:
            paths.append(part)
    return paths


def choose_column(df, preferred, candidates):
    if preferred in df.columns:
        return preferred
    lowered = {str(col).strip().lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    return None


def build_annotation_map(csv_paths, csv_cell_col, csv_label_col):
    mappings = []
    for csv_path in csv_paths:
        if not os.path.isfile(csv_path):
            sys.exit(f"annotation csv not found: {csv_path}")
        df = pd.read_csv(csv_path)
        if df.shape[0] == 0:
            continue
        cell_col = choose_column(df, csv_cell_col, ["barcode", "cell_id", "cellid", "cell_barcode"])
        label_col = choose_column(df, csv_label_col, ["grouped_annotation", "celltype", "annotation", "cluster_id", "group"])
        if cell_col is None or label_col is None:
            sys.exit(f"required columns not found in {csv_path}. got columns: {','.join(df.columns.astype(str))}")
        sub = df[[cell_col, label_col]].copy()
        sub.columns = ["cell_key", "annotation"]
        sub["cell_key"] = sub["cell_key"].astype(str)
        sub["annotation"] = sub["annotation"].astype(str)
        sub = sub[sub["cell_key"].str.len() > 0]
        mappings.append(sub)
    if len(mappings) == 0:
        return {}
    merged = pd.concat(mappings, axis=0, ignore_index=True)
    merged = merged.drop_duplicates(subset=["cell_key"], keep="last")
    return dict(zip(merged["cell_key"], merged["annotation"]))


def merge_reannotation_to_base(
    base_zarr,
    output_dir,
    annotation_csv,
    csv_cell_col,
    csv_label_col,
    input_cell_col,
    target_col,
    fallback_col
):
    csv_paths = parse_csv_inputs(annotation_csv)
    if len(csv_paths) == 0:
        sys.exit("merge_by=reannotation requires --annotation_csv")

    sdata = spd.read_zarr(base_zarr)
    table_keys = list(sdata.tables.keys())
    if len(table_keys) == 0:
        sys.exit("no tables found in base zarr")
    table_key = table_keys[0]
    table = sdata.tables[table_key]

    mapping = build_annotation_map(csv_paths, csv_cell_col, csv_label_col)
    if len(mapping) == 0:
        sys.exit("no valid annotation rows found in csv inputs")

    if input_cell_col in table.obs.columns:
        cell_series = table.obs[input_cell_col].astype(str)
    else:
        cell_series = table.obs.index.astype(str)

    if fallback_col in table.obs.columns:
        merged_labels = table.obs[fallback_col].astype(str).copy()
    elif target_col in table.obs.columns:
        merged_labels = table.obs[target_col].astype(str).copy()
    else:
        merged_labels = pd.Series(["Unknown"] * table.n_obs, index=table.obs.index)

    mapped = cell_series.map(mapping)
    mask = mapped.notna()
    merged_labels.loc[mask] = mapped.loc[mask].astype(str)
    table.obs[target_col] = merged_labels.astype("category")

    table.obs["cell_id"] = table.obs.index.astype(str)
    table.uns["spatialdata_attrs"]["instance_key"] = "cell_id"
    sdata.tables[table_key] = table
    os.makedirs(output_dir, exist_ok=True)
    sdata.write(os.path.join(output_dir, "concatenated_sdata.zarr"), overwrite=True)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge')
    parser.add_argument('--INPUT', nargs='+', required=True, 
                       help='INPUT zarr')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='output_dir')
    parser.add_argument('--merge_by', type=str, 
                       required=True)
    parser.add_argument('--reordering', type=str, default='False',
                       help='重排聚类标签')
    parser.add_argument('--re_sample', type=str, default='False',
                       help='样本标识')
    parser.add_argument('--cluster_key', type=str, default='leiden',
                       help='聚类列名')
    parser.add_argument('--annotation_csv', type=str, default='',
                       help='csv path, directory, or comma-separated csv paths')
    parser.add_argument('--csv_cell_col', type=str, default='Barcode',
                       help='cell id column in csv')
    parser.add_argument('--csv_label_col', type=str, default='Grouped_Annotation',
                       help='annotation column in csv')
    parser.add_argument('--input_cell_col', type=str, default='cell_id',
                       help='cell id column in base zarr table obs')
    parser.add_argument('--target_col', type=str, default='celltype',
                       help='target column to write merged annotation')
    parser.add_argument('--fallback_col', type=str, default='celltype',
                       help='fallback column for cells not in csv')
    args = parser.parse_args()
    print("starting .......................")
    print(args.INPUT)
    reordering = args.reordering.lower() == 'true'
    re_sample = args.re_sample.lower() == 'true'
    if args.merge_by == "sample":
        sdatas = merge_by_samples(args.INPUT, re_sample)
        concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
        concatenated_sdata.write(os.path.join(args.output_dir,"concatenated_sdata.zarr"), overwrite=True)
        del concatenated_sdata, sdatas
    elif args.merge_by == "reannotation":
        merge_reannotation_to_base(
            base_zarr=args.INPUT[0],
            output_dir=args.output_dir,
            annotation_csv=args.annotation_csv,
            csv_cell_col=args.csv_cell_col,
            csv_label_col=args.csv_label_col,
            input_cell_col=args.input_cell_col,
            target_col=args.target_col,
            fallback_col=args.fallback_col
        )
    else:
        sdatas = merge_by_clusters(args.INPUT, reordering, args.cluster_key)
        concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
        concatenated_sdata.write(os.path.join(args.output_dir,"concatenated_sdata.zarr"), overwrite=True)
        del concatenated_sdata, sdatas
    gc.collect()











