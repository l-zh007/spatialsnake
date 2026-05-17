import os
import warnings
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import spatialdata as spd
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv

warnings.filterwarnings("ignore")


def load_spatialdata(input_path):
    if not input_path.endswith(".zarr"):
        raise ValueError("reclustering input must be a .zarr path")
    sdata = spd.read_zarr(input_path)
    table_keys = list(sdata.tables.keys())
    if len(table_keys) == 0:
        raise ValueError("no tables found in input zarr")
    table_key = table_keys[0]
    adata = sdata[table_key].copy()
    return sdata, table_key, adata


def run_reclustering(adata, resolution, n_top_genes, neighbors, n_pcs, cluster_key):
    #sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
    #if "highly_variable" in adata.var.columns and adata.var["highly_variable"].sum() > 0:
        #adata = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(adata, svd_solver="arpack")
    n_pcs_use = min(int(n_pcs), int(adata.obsm["X_pca"].shape[1]))
    sc.pp.neighbors(adata, n_neighbors=int(neighbors), n_pcs=n_pcs_use)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=float(resolution), key_added=cluster_key, flavor="igraph", random_state=0)
    return adata


def export_marker_genes(adata, cluster_key, method, min_pct, logfc_threshold, output_csv):
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=method,
            min_in_group_fraction=float(min_pct),
            log2fc_min=float(logfc_threshold)
        )
    except TypeError:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method)
    groups = adata.obs[cluster_key].astype(str).unique().tolist()
    frames = []
    for group_name in groups:
        df = sc.get.rank_genes_groups_df(adata, group=group_name)
        if "pct_nz_group" in df.columns:
            df = df[df["pct_nz_group"] >= float(min_pct)]
        if "logfoldchanges" in df.columns:
            df = df[df["logfoldchanges"] >= float(logfc_threshold)]
            df = df.rename(columns={"logfoldchanges": "avg_log2FC"})
        else:
            df["avg_log2FC"] = np.nan
        df["cluster"] = str(group_name)
        if "pvals" in df.columns:
            df = df.rename(columns={"pvals": "p_val"})
        elif "pval" not in df.columns:
            df["p_val"] = np.nan
        col_order = ["cluster", "names", "p_val", "avg_log2FC", "pvals_adj", "scores", "pct_nz_group", "pct_nz_reference"]
        available_cols = [c for c in col_order if c in df.columns]
        df = df[available_cols].copy()
        if "names" in df.columns:
            df = df.rename(columns={"names": "gene"})
        frames.append(df)
    marker_df = pd.concat(frames, axis=0, ignore_index=True) if len(frames) > 0 else pd.DataFrame(columns=["cluster", "gene", "p_val", "avg_log2FC"])
    marker_df.to_csv(output_csv, index=False)


def save_umap_plot(adata, cluster_key, output_png):
    sc.pl.umap(adata, color=[cluster_key], show=False)
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()


def save_spatial_plot(sdata, table_key, adata, cluster_key, output_png):
    sdata[table_key] = adata
    try:
        coord_systems = sorted(list(sdata.coordinate_systems))
        coord = coord_systems[0] if len(coord_systems) > 0 else None
        images = list(sdata.images.keys())
        if len(images) > 0 and coord is not None:
            ax = plt.subplots(1, 1, figsize=(8, 8))[1]
            sdata.pl.render_images().pl.render_shapes(color=cluster_key).pl.show(ax=ax, coordinate_systems=coord)
            plt.savefig(output_png, dpi=300, bbox_inches="tight")
            plt.close()
            return
    except Exception:
        pass
    if "spatial" in adata.obsm_keys():
        coords = adata.obsm["spatial"]
        codes = adata.obs[cluster_key].astype("category").cat.codes
        plt.figure(figsize=(8, 8))
        plt.scatter(coords[:, 0], coords[:, 1], c=codes, s=5, cmap="tab20", linewidths=0)
        plt.gca().invert_yaxis()
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        sc.pl.umap(adata, color=[cluster_key], show=False)
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="spatial reclustering")
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--sample_id", type=str, required=True)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--marker_method", type=str, default="wilcoxon")
    parser.add_argument("--min_pct", type=float, default=0.1)
    parser.add_argument("--logfc_threshold", type=float, default=0.25)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    cluster_key = "recluster"
    sdata, table_key, adata = load_spatialdata(args.input)
    print(adata)
    print(adata.X)
    print(adata.X.shape)
    adata = run_reclustering(adata, args.resolution, args.n_top_genes, args.neighbors, args.n_pcs, cluster_key)
    umap_png = os.path.join(args.output_dir, "umap_recluster.png")
    spatial_png = os.path.join(args.output_dir, "spatial_clusters.png")
    marker_csv = os.path.join(args.output_dir, "marker_genes.csv")
    output_zarr = os.path.join(args.output_dir, f"{args.sample_id}.zarr")

    save_umap_plot(adata, cluster_key, umap_png)
    save_spatial_plot(sdata, table_key, adata, cluster_key, spatial_png)
    export_marker_genes(adata, cluster_key, args.marker_method, args.min_pct, args.logfc_threshold, marker_csv)
    exported_files = export_cluster_csv(
        adata,
        data_type="visium",
        dir_path=args.output_dir,
        cell_id_col="cell_id",
        info_col=cluster_key,
        sample_col="sample",
        sample_id=args.sample_id
    )
    if len(exported_files) > 0:
        os.replace(exported_files[0], os.path.join(args.output_dir, "cluster_assignments.csv"))

    sdata[table_key] = adata
    sdata.write(output_zarr, overwrite=True)


if __name__ == "__main__":
    main()
