import os
import sys
import warnings
import argparse


def _thread_count_from_argv(default=8):
    """Read --threads early so BLAS is configured before NumPy/Scanpy import."""
    argv = sys.argv[1:]
    for index, token in enumerate(argv):
        raw_value = None
        if token == "--threads" and index + 1 < len(argv):
            raw_value = argv[index + 1]
        elif token.startswith("--threads="):
            raw_value = token.split("=", 1)[1]
        if raw_value is not None:
            try:
                value = int(raw_value)
            except (TypeError, ValueError):
                return default
            return value if value >= 1 else default
    return default


def configure_thread_environment(threads):
    threads = int(threads)
    if threads < 1:
        raise ValueError("threads must be >= 1")
    value = str(threads)
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = value
    return threads


configure_thread_environment(_thread_count_from_argv())

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import spatialdata as spd
import spatialdata_plot
from spatialdata.transformations import get_transformation
from spatialsnake.workflow.function.export_cluster_csv import export_cluster_csv
from spatialsnake.workflow.function.stereoseq_selection import resolve_region_image
from spatialsnake.workflow.function.logging_utils import setup_logger, log_step

warnings.filterwarnings("ignore")
logger = setup_logger("reclustering")


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
    if int(n_top_genes) < 0 or int(n_top_genes) == 1 or min(adata.shape) < 3:
        raise ValueError("Reclustering requires at least 3 observations/genes; n_top_genes must be 0 (all genes) or >= 2")
    use_hvg = int(n_top_genes) > 0
    if use_hvg:
        effective_top = min(int(n_top_genes), adata.n_vars)
        if effective_top == adata.n_vars:
            adata.var['highly_variable'] = True
        else:
            sc.pp.highly_variable_genes(adata, n_top_genes=effective_top, flavor="seurat")
        if int(adata.var['highly_variable'].sum()) < 2:
            raise ValueError("Fewer than two subset highly variable genes remain; increase n_top_genes")
        logger.info(f"Reclustering uses {int(adata.var['highly_variable'].sum())} subset highly variable genes")
    else:
        logger.info(f"HVG selection disabled; reclustering PCA uses all {adata.n_vars} genes")
    sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=use_hvg)
    n_pcs_use = min(int(n_pcs), int(adata.obsm["X_pca"].shape[1]))
    sc.pp.neighbors(adata, n_neighbors=int(neighbors), n_pcs=n_pcs_use)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=float(resolution), key_added=cluster_key, flavor="igraph", random_state=0)
    return adata


def export_marker_genes(adata, cluster_key, method, min_pct, logfc_threshold, output_csv):
    sc.tl.rank_genes_groups(
        adata, groupby=cluster_key, method=method, use_raw=False, pts=True
    )
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
        attrs = adata.uns.get("spatialdata_attrs", {})
        regions = attrs.get("region", [])
        regions = [regions] if isinstance(regions, str) else list(regions or [])
        regions = [region for region in regions if region in sdata.shapes or region in sdata.points]
        images = [name for name in sdata.images if "hires" in name.lower()] or list(sdata.images)
        systems = sorted(sdata.coordinate_systems)
        panels = []
        for region in regions:
            element = sdata.shapes[region] if region in sdata.shapes else sdata.points[region]
            candidates = [system for system in systems if system in get_transformation(element, get_all=True)]
            if not candidates:
                continue
            system = next((system for system in candidates if "hires" in system.lower()), candidates[0])
            image = resolve_region_image(sdata, adata, region, system)
            panels.append((region, system, image))

        if panels:
            rows = 2 if any(image for _, _, image in panels) else 1
            _, axes = plt.subplots(rows, len(panels), figsize=(10 * len(panels), 13 if rows == 2 else 8), squeeze=False)
            groups = adata.obs[cluster_key].astype("category").cat.categories.tolist()
            instance_key = attrs.get("instance_key")
            region_key = attrs.get("region_key")
            for i, (region, system, image) in enumerate(panels):
                ax = axes[-1, i]
                if image:
                    sdata.pl.render_images(image).pl.show(ax=axes[0, i], coordinate_systems=system, title="image")
                if region in sdata.shapes:
                    sdata.pl.render_shapes(
                        region, color="#9E9E9E", fill_alpha=1, outline_alpha=0,
                    ).pl.render_shapes(
                        region, color=cluster_key, groups=groups, table_name=table_key,
                        na_color=None, fill_alpha=1, outline_alpha=0,
                    ).pl.show(ax=ax, coordinate_systems=system, title=str(region))
                else:
                    sdata.pl.render_points(region, color="#9E9E9E", size=2, alpha=1).pl.show(
                        ax=ax, coordinate_systems=system, title=str(region),
                    )
                    region_obs = adata.obs
                    if region_key in region_obs:
                        region_obs = region_obs[region_obs[region_key].astype(str) == str(region)]
                    label_ids = region_obs[instance_key].astype(str) if instance_key in region_obs else region_obs.index.astype(str)
                    labels_by_id = pd.Series(region_obs[cluster_key].to_numpy(), index=label_ids)
                    points = spd.transform(sdata.points[region], to_coordinate_system=system).compute()
                    points.index = points.index.astype(str)
                    labels = labels_by_id.reindex(points.index).astype("category").cat.remove_unused_categories()
                    colors = plt.get_cmap("tab20", max(1, len(labels.cat.categories)))
                    for j, label in enumerate(labels.cat.categories):
                        selected = labels == label
                        ax.scatter(points.loc[selected, "x"], points.loc[selected, "y"], s=3, color=colors(j), label=label)
                    if len(labels.cat.categories):
                        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
            plt.tight_layout()
            plt.savefig(output_png, dpi=300, bbox_inches="tight")
            plt.close()
            return
    except Exception as exc:
        logger.warning("Spatial context plotting failed; using coordinate fallback: %s", exc)
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
    parser.add_argument("--n_top_genes", type=int, default=0)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--marker_method", type=str, default="wilcoxon")
    parser.add_argument("--min_pct", type=float, default=0.1)
    parser.add_argument("--logfc_threshold", type=float, default=0.25)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be >= 1")
    configure_thread_environment(args.threads)
    sc.settings.n_jobs = args.threads

    os.makedirs(args.output_dir, exist_ok=True)
    cluster_key = "recluster"
    logger.info(f"Using {args.threads} threads for sample {args.sample_id}")
    log_step(logger, 1, 5, "loading SpatialData input")
    sdata, table_key, adata = load_spatialdata(args.input)
    logger.info(f"Loaded {adata.n_obs} observations and {adata.n_vars} genes")
    log_step(logger, 2, 5, "running PCA, neighbors, UMAP and Leiden reclustering")
    adata = run_reclustering(adata, args.resolution, args.n_top_genes, args.neighbors, args.n_pcs, cluster_key)
    umap_png = os.path.join(args.output_dir, "umap_recluster.png")
    spatial_png = os.path.join(args.output_dir, "spatial_clusters.png")
    marker_csv = os.path.join(args.output_dir, "marker_genes.csv")
    output_zarr = os.path.join(args.output_dir, f"{args.sample_id}.zarr")

    log_step(logger, 3, 5, "saving UMAP and spatial plots")
    save_umap_plot(adata, cluster_key, umap_png)
    save_spatial_plot(sdata, table_key, adata, cluster_key, spatial_png)
    log_step(logger, 4, 5, "exporting marker genes and cluster assignments")
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
    log_step(logger, 5, 5, f"saving reclustered data to {output_zarr}")
    sdata.write(output_zarr, overwrite=True)
    logger.info("Reclustering module completed")


if __name__ == "__main__":
    main()
