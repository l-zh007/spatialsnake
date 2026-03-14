import os
import pandas as pd
import spatialdata as spd
import scanpy as sc

def export_cluster_csv(data_input, data_type, dir_path, cell_id_col="cell_id", info_col="clusters", sample_col="sample", sample_id=None):
    os.makedirs(dir_path, exist_ok=True)
    adata = None
    if isinstance(data_input, str):
        if data_input.endswith(".zarr"):
            sdata = spd.read_zarr(data_input)
            table_names = list(sdata.tables.keys())
            if len(table_names) == 0:
                raise ValueError("no tables found in zarr input")
            adata = sdata[table_names[0]]
        elif data_input.endswith(".h5ad"):
            adata = sc.read_h5ad(data_input)
        else:
            raise ValueError("data_input must be .zarr or .h5ad path")
    else:
        if hasattr(data_input, "tables"):
            table_names = list(data_input.tables.keys())
            if len(table_names) == 0:
                raise ValueError("no tables found in spatialdata input")
            adata = data_input[table_names[0]]
        elif hasattr(data_input, "obs"):
            adata = data_input
        else:
            raise ValueError("data_input must be SpatialData, AnnData, or path")

    if info_col not in adata.obs.columns:
        if "grouped_clusters" in adata.obs.columns:
            info_col = "grouped_clusters"
        elif "clusters" in adata.obs.columns:
            info_col = "clusters"
        else:
            raise ValueError(f"{info_col} not found in obs")

    effective_sample_col = sample_col if sample_col in adata.obs.columns else None
    if effective_sample_col is None and "region" in adata.obs.columns:
        effective_sample_col = "region"

    if effective_sample_col:
        sample_values = [v for v in adata.obs[effective_sample_col].unique() if str(v) != "nan"]
        if len(sample_values) == 0:
            sample_values = [sample_id or "sample"]
    else:
        sample_values = [sample_id or "sample"]

    data_type_value = str(data_type).lower() if data_type is not None else ""
    use_xenium = "xenium" in data_type_value

    output_files = []
    for sample_name in sample_values:
        sample_name = str(sample_name)
        if effective_sample_col:
            adata_sample = adata[adata.obs[effective_sample_col] == sample_name].copy()
        else:
            adata_sample = adata.copy()

        cell_id_series = adata_sample.obs[cell_id_col] if cell_id_col in adata_sample.obs.columns else adata_sample.obs_names
        cell_id_series = cell_id_series.astype(str)
        group_series = adata_sample.obs[info_col].astype(str)

        if use_xenium:
            df_output = pd.DataFrame({
                "cell_id": cell_id_series,
                "group": group_series
            })
            output_filename = os.path.join(dir_path, f"{sample_name}_cell_groups.csv")
        else:
            barcode_series = cell_id_series
            if barcode_series.str.contains("cellid_").any():
                barcode_series = "cellid_" + barcode_series.str.split("cellid_").str[-1]
            df_output = pd.DataFrame({
                "Barcode": barcode_series,
                "Grouped_Annotation": group_series
            })
            output_filename = os.path.join(dir_path, f"{sample_name}_cell_clusters.csv")

        df_output.to_csv(output_filename, index=False)
        output_files.append(output_filename)

    return output_files
