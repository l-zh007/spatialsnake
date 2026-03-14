import os
import json
import argparse
import scanpy as sc
import spatialdata as spd
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv


def resolve_annotation_map(anno_data, sample_id):
    if sample_id in anno_data:
        return anno_data[sample_id]
    first_key = next(iter(anno_data.keys()))
    return anno_data[first_key]


def annotate_adata(adata, annotation_map):
    cluster_key = "recluster" if "recluster" in adata.obs.columns else "clusters"
    if cluster_key not in adata.obs.columns:
        raise ValueError("no recluster or clusters column found in table obs")
    mapped = adata.obs[cluster_key].astype(str).map(annotation_map)
    adata.obs["celltype"] = mapped.fillna("Unknown").astype("category")
    return adata


def main():
    parser = argparse.ArgumentParser(description="reannotation")
    parser.add_argument("--input_dir", type=str, required=True)
    parser.add_argument("--sample_id", type=str, required=True)
    parser.add_argument("--output_zarr_path", type=str, required=True)
    parser.add_argument("--output_csv", type=str, required=True)
    parser.add_argument("--type", type=str, required=True)
    parser.add_argument("--anno_data", type=str, required=True)
    args = parser.parse_args()

    anno_data = json.loads(args.anno_data)
    annotation_map = resolve_annotation_map(anno_data, args.sample_id)
    output_dir = os.path.dirname(args.output_csv)
    os.makedirs(output_dir, exist_ok=True)

    if args.input_dir.endswith(".h5ad"):
        adata = sc.read_h5ad(args.input_dir)
        adata = annotate_adata(adata, annotation_map)
        exported_files = export_cluster_csv(
            adata,
            args.type,
            output_dir,
            cell_id_col="cell_id",
            info_col="celltype",
            sample_col="region",
            sample_id=args.sample_id
        )
        adata.write(args.output_zarr_path)
    else:
        sdata = spd.read_zarr(args.input_dir)
        table_keys = list(sdata.tables.keys())
        if len(table_keys) == 0:
            raise ValueError("no tables found in input zarr")
        table_key = table_keys[0]
        adata = sdata[table_key].copy()
        adata = annotate_adata(adata, annotation_map)
        sdata[table_key] = adata
        exported_files = export_cluster_csv(
            sdata,
            args.type,
            output_dir,
            cell_id_col="cell_id",
            info_col="celltype",
            sample_col="region",
            sample_id=args.sample_id
        )
        sdata.write(args.output_zarr_path, overwrite=True)

    if len(exported_files) > 0:
        os.replace(exported_files[0], args.output_csv)


if __name__ == "__main__":
    main()
