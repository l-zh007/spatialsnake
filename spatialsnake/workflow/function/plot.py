import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import os
import re
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import gc
from matplotlib.colors import TwoSlopeNorm
from scipy import sparse
from spatialsnake.workflow.function.logging_utils import setup_logger

logger = setup_logger("plot")

def cluster_proportion(adata, sample_col='region', 
                      cluster_col='clusters',
                      figsize=(12, 8),
                      palette='tab20',
                      sort_samples=None,
                      sort_clusters=None,
                      title=None,
                      save_path=None,
                      dpi=300):
    """
    绘制堆叠条形图，展示每个样本中各cluster的比例分布
    
    Parameters:
    adata: AnnData对象 - 包含空间转录组数据
    sample_col: str - 样本列名
    cluster_col: str - cluster列名  
    figsize: tuple - 图表大小
    palette: str or list - 颜色调色板
    sort_samples: list - 样本排序顺序
    sort_clusters: list - cluster排序顺序
    title: str - 图表标题
    save_path: str - 保存路径
    dpi: int - 图片分辨率
    """
    
    # 数据验证
    if sample_col not in adata.obs.columns:
        raise ValueError(f"在obs中没有找到'{sample_col}'列")
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"在obs中没有找到'{cluster_col}'列")
    
    # 计算每个样本中各cluster的数量
    cluster_counts = pd.crosstab(adata.obs[sample_col], adata.obs[cluster_col])
    
    # 计算比例
    cluster_props = cluster_counts.div(cluster_counts.sum(axis=1), axis=0) * 100
    
    # 排序
    if sort_samples is not None:
        cluster_props = cluster_props.reindex(sort_samples)
    if sort_clusters is not None:
        cluster_props = cluster_props[sort_clusters]
    
    # 创建图表
    fig, ax = plt.subplots(figsize=figsize)
    
    # 绘制堆叠条形图
    cluster_props.plot(kind='bar', stacked=True, ax=ax, 
                      colormap=palette, width=0.7)
    
    # 美化图表
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Proportion (%)', fontsize=12, fontweight='bold')
    
    if title is None:
        title = f'Cluster Proportion Distribution by Sample'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # 设置y轴范围
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.set_yticklabels([f'{x}%' for x in range(0, 101, 20)], fontsize=10)
    
    # 美化x轴标签
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    
    # 美化图例
    ax.legend(title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left', 
             fontsize=9, title_fontsize=10)
    
    # 添加网格
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {save_path}")
    
    return fig, ax
  

def plot_celltype_spatial_enrichment(
    adata,
    celltype_col="celltype",
    cluster_col="spatial_cluster",
    save_path=None,
    title=None,
    min_count=1,
    annotate=True,
    cmap="PuOr",
    dpi=320,
):
    """
    绘制 celltype 在各 spatial cluster 中的富集热图。

    热图数值使用 log2(observed / expected)，其中 expected 基于列联表的独立性假设计算。
    横坐标为 celltype，纵坐标为 spatial_cluster。
    """

    if celltype_col not in adata.obs.columns:
        raise ValueError(f"在obs中没有找到'{celltype_col}'列")
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"在obs中没有找到'{cluster_col}'列")

    plot_df = adata.obs[[celltype_col, cluster_col]].copy()
    plot_df = plot_df.dropna()
    if plot_df.shape[0] == 0:
        raise ValueError("用于绘图的数据为空，无法计算富集热图")

    contingency = pd.crosstab(plot_df[cluster_col].astype(str), plot_df[celltype_col].astype(str))
    contingency = contingency.loc[contingency.sum(axis=1) >= min_count, contingency.sum(axis=0) >= min_count]
    if contingency.shape[0] == 0 or contingency.shape[1] == 0:
        raise ValueError("过滤低频 cluster/celltype 后无有效数据可绘图")

    total = contingency.values.sum()
    expected = np.outer(contingency.sum(axis=1), contingency.sum(axis=0)) / total
    observed = contingency.to_numpy(dtype=float)
    enrichment = np.log2((observed + 1e-6) / (expected + 1e-6))
    enrichment_df = pd.DataFrame(
        enrichment,
        index=contingency.index.astype(str),
        columns=contingency.columns.astype(str),
    )

    ordered_clusters = (
        enrichment_df.abs().max(axis=1).sort_values(ascending=False).index.tolist()
    )
    ordered_celltypes = (
        enrichment_df.abs().max(axis=0).sort_values(ascending=False).index.tolist()
    )
    enrichment_df = enrichment_df.loc[ordered_clusters, ordered_celltypes]

    fig_width = max(7.5, 0.7 * enrichment_df.shape[1] + 2.6)
    fig_height = max(6.5, 0.52 * enrichment_df.shape[0] + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

    vmax = float(np.nanmax(np.abs(enrichment_df.to_numpy())))
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(
        enrichment_df.to_numpy(),
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(enrichment_df.shape[1]))
    ax.set_xticklabels(enrichment_df.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(np.arange(enrichment_df.shape[0]))
    ax.set_yticklabels(enrichment_df.index, fontsize=10)
    ax.set_xlabel("Cell type", fontsize=12, labelpad=10)
    ax.set_ylabel("Spatial cluster", fontsize=12, labelpad=10)
    ax.set_facecolor("white")

    if title is None:
        title = "Cell type enrichment across spatial clusters"
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)

    if annotate and enrichment_df.shape[0] * enrichment_df.shape[1] <= 180:
        for row_idx in range(enrichment_df.shape[0]):
            for col_idx in range(enrichment_df.shape[1]):
                value = enrichment_df.iat[row_idx, col_idx]
                ax.text(
                    col_idx,
                    row_idx,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color="white" if abs(value) > vmax * 0.45 else "black",
                    fontweight="bold",
                )

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log2(observed / expected)", fontsize=10, labelpad=8)
    cbar.ax.tick_params(labelsize=9, length=2)
    cbar.outline.set_linewidth(0.5)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xticks(np.arange(-0.5, enrichment_df.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, enrichment_df.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#f0f0f0", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()

    if save_path is not None:
        save_path = str(save_path)
        if not save_path.lower().endswith(".png"):
            save_path = f"{save_path}.png"
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", pad_inches=0.08, facecolor="white")
        logger.info(f"Plot saved to {save_path}")
    return fig, ax, enrichment_df

def plot_auc_heatmap_scanpy(auc_mtx,adata, groupby, top_genes=None, figsize=(10, 8),outputs="results/scenic_results/heatmap.png"):
    """
    使用Scanpy绘制AUCell热图
    """
    # 创建临时AnnData对象用于热图
    temp_adata = sc.AnnData(X=auc_mtx)  # 转置以匹配scanpy的格式（基因×细胞）
    temp_adata.obs = adata.obs.copy()
    temp_adata.var_names = auc_mtx.columns  # TF名称

    if top_genes is not None:
        mean_auc = auc_mtx.mean().sort_values(ascending=False)
        selected_tfs = mean_auc.head(top_genes).index.tolist()
        temp_adata = temp_adata[:, selected_tfs]
    
    # 绘制热图
    sc.pl.heatmap(
        temp_adata,
        var_names=temp_adata.var_names.tolist(),
        groupby=groupby,
        figsize=figsize,
        cmap='Reds',  # 使用红色系，适合AUC值
        dendrogram=False,
        standard_scale='var',  # 按TF标准化
        swap_axes=True,  # 交换轴，使TF在y轴
        show_gene_labels=True,
        show=False)
    plt.savefig(outputs, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
def plot_cell2location_dotplot(
    adata,
    unsupervised_cluster_col,
    save_path,
    sample_col=None,
    abundance_key="q05_cell_abundance_w_sf",
    max_cell_types=30,
    enrichment_clip=2.5,
    source_data_path=None,
    size_legend_title=None,
    plot_title=None,
    x_axis_label="Unsupervised spatial cluster",
    method_label="Cell2location",
):
    """Plot sample-balanced regional enrichment of continuous cell2location abundance.

    Bubble area encodes the mean spot-level abundance within a spatial cluster.
    Colour encodes the log2 ratio between that cluster mean and the tissue-wide
    mean for the same cell type. For multiple samples, both quantities are first
    calculated within each sample and then averaged with equal sample weight.
    """

    def _natural_key(value):
        return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(value))]

    def _cluster_order(series):
        if isinstance(series.dtype, pd.CategoricalDtype):
            observed = set(series.dropna().astype(str))
            return [str(value) for value in series.cat.categories if str(value) in observed]
        return sorted(series.dropna().astype(str).unique().tolist(), key=_natural_key)

    if unsupervised_cluster_col not in adata.obs.columns:
        raise KeyError(f"'{unsupervised_cluster_col}' was not found in adata.obs")
    if abundance_key not in {"q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"}:
        raise ValueError("abundance_key must be q05_cell_abundance_w_sf or means_cell_abundance_w_sf")
    if abundance_key not in adata.obsm:
        raise KeyError(f"Expected cell2location abundance matrix missing from obsm: {abundance_key}")
    if max_cell_types is not None and int(max_cell_types) < 0:
        raise ValueError("max_cell_types must be >= 0; use 0 to display all cell types")
    if not np.isfinite(enrichment_clip) or enrichment_clip <= 0:
        raise ValueError("enrichment_clip must be a positive finite number")

    factor_names = adata.uns.get("mod", {}).get("factor_names")
    if factor_names is None:
        raise KeyError("Expected cell2location cell-type names in uns['mod']['factor_names']")
    factor_names = pd.Index([str(value) for value in factor_names])
    if factor_names.empty or factor_names.has_duplicates or (factor_names.str.strip() == "").any():
        raise ValueError("cell2location factor_names must be non-empty and unique")

    abundance = adata.obsm[abundance_key]
    if sparse.issparse(abundance):
        abundance = abundance.toarray()
    if isinstance(abundance, pd.DataFrame):
        abundance_df = abundance.copy()
        abundance_df.index = abundance_df.index.astype(str)
        target_index = pd.Index(adata.obs_names.astype(str))
        if abundance_df.index.is_unique and set(abundance_df.index) == set(target_index):
            abundance_df = abundance_df.loc[target_index]
        elif len(abundance_df) == len(target_index) and isinstance(abundance_df.index, pd.RangeIndex):
            abundance_df.index = target_index
        elif not abundance_df.index.equals(target_index):
            raise ValueError(f"{abundance_key} rows do not align with adata.obs_names")
    else:
        abundance = np.asarray(abundance)
        if abundance.ndim != 2 or abundance.shape[0] != adata.n_obs:
            raise ValueError(
                f"Invalid {abundance_key} shape {abundance.shape}; expected ({adata.n_obs}, n_cell_types)"
            )
        abundance_df = pd.DataFrame(abundance, index=adata.obs_names.astype(str))
    if abundance_df.shape[1] != len(factor_names):
        raise ValueError(
            f"{abundance_key} has {abundance_df.shape[1]} columns but factor_names has {len(factor_names)}"
        )
    abundance_df.columns = factor_names
    abundance_df = abundance_df.astype(float)
    abundance_values = abundance_df.to_numpy()
    if not np.all(np.isfinite(abundance_values)) or np.any(abundance_values < 0):
        raise ValueError(f"{abundance_key} must contain finite, non-negative values")

    cluster_series = adata.obs[unsupervised_cluster_col].copy()
    cluster_order = _cluster_order(cluster_series)
    cluster_text = cluster_series.astype("string").str.strip()
    valid_mask = cluster_text.notna() & cluster_text.ne("")

    if sample_col is not None:
        if sample_col not in adata.obs.columns:
            raise KeyError(f"Sample column '{sample_col}' was not found in adata.obs")
        sample_series = adata.obs[sample_col].astype("string").str.strip()
        valid_mask &= sample_series.notna() & sample_series.ne("")
    else:
        sample_series = pd.Series("sample", index=adata.obs_names, dtype="string")

    if not valid_mask.any():
        raise ValueError("No spots with valid cluster and sample labels were available for the dotplot")
    abundance_df = abundance_df.loc[valid_mask.to_numpy()].copy()
    cluster_text = cluster_text.loc[valid_mask].astype(str)
    sample_series = sample_series.loc[valid_mask].astype(str)
    cluster_order = [value for value in cluster_order if value in set(cluster_text)]

    metadata = pd.DataFrame(
        {"sample": sample_series.to_numpy(), "cluster": cluster_text.to_numpy()},
        index=abundance_df.index,
    )
    grouped_input = abundance_df.copy()
    grouped_input[["__sample", "__cluster"]] = metadata[["sample", "cluster"]]
    sample_cluster_mean = grouped_input.groupby(
        ["__sample", "__cluster"], observed=True, sort=False
    )[factor_names].mean()
    sample_tissue_mean = grouped_input.groupby("__sample", observed=True, sort=False)[factor_names].mean()

    global_mean = sample_tissue_mean.mean(axis=0)
    keep_cell_types = global_mean[global_mean > 0].index
    if len(keep_cell_types) == 0:
        raise ValueError(f"All cell types in {abundance_key} have zero abundance")
    sample_cluster_mean = sample_cluster_mean.loc[:, keep_cell_types]
    sample_tissue_mean = sample_tissue_mean.loc[:, keep_cell_types]
    global_mean = global_mean.loc[keep_cell_types]

    sample_ids = sample_cluster_mean.index.get_level_values("__sample")
    denominator = sample_tissue_mean.loc[sample_ids].to_numpy(dtype=float)
    numerator = sample_cluster_mean.to_numpy(dtype=float)
    pseudocount = np.where(denominator > 0, denominator * 1e-6, np.nan)
    sample_enrichment = pd.DataFrame(
        np.log2((numerator + pseudocount) / (denominator + pseudocount)),
        index=sample_cluster_mean.index,
        columns=sample_cluster_mean.columns,
    )

    mean_abundance = sample_cluster_mean.groupby(level="__cluster", sort=False).mean().T
    log2_enrichment = sample_enrichment.groupby(level="__cluster", sort=False).mean().T
    mean_abundance = mean_abundance.reindex(columns=cluster_order)
    log2_enrichment = log2_enrichment.reindex(columns=cluster_order)

    cluster_spot_counts = metadata.groupby("cluster", sort=False).size().reindex(cluster_order)
    cluster_sample_counts = metadata.groupby("cluster", sort=False)["sample"].nunique().reindex(cluster_order)
    total_samples = metadata["sample"].nunique()
    if total_samples > 1:
        low_replication = cluster_sample_counts[cluster_sample_counts < 2].index.tolist()
        if low_replication:
            logger.warning(
                "Regional enrichment is descriptive: %d cluster(s) occur in only one sample. Examples: %s",
                len(low_replication),
                low_replication[:10],
            )
    ranked_cell_types = global_mean.sort_values(ascending=False, kind="stable").index.tolist()
    max_cell_types = int(max_cell_types) if max_cell_types is not None else 0
    plotted_cell_types = ranked_cell_types if max_cell_types == 0 else ranked_cell_types[:max_cell_types]
    omitted_cell_types = ranked_cell_types[len(plotted_cell_types):]
    if omitted_cell_types:
        logger.warning(
            "Dotplot displays the top %d of %d cell types by sample-balanced tissue mean; "
            "all values remain in the source TSV",
            len(plotted_cell_types),
            len(ranked_cell_types),
        )

    peak_cluster = log2_enrichment.loc[plotted_cell_types].idxmax(axis=1)
    cluster_positions = {cluster: index for index, cluster in enumerate(cluster_order)}
    cell_type_order = sorted(
        plotted_cell_types,
        key=lambda cell_type: (
            cluster_positions.get(peak_cluster.get(cell_type), len(cluster_order)),
            -float(global_mean[cell_type]),
            _natural_key(cell_type),
        ),
    )
    plot_abundance = mean_abundance.loc[cell_type_order, cluster_order]
    plot_enrichment = log2_enrichment.loc[cell_type_order, cluster_order]

    source_data = (
        mean_abundance.rename_axis(index="cell_type", columns="cluster")
        .reset_index()
        .melt(id_vars="cell_type", var_name="cluster", value_name="mean_abundance")
        .merge(
            log2_enrichment.rename_axis(index="cell_type", columns="cluster")
            .reset_index()
            .melt(id_vars="cell_type", var_name="cluster", value_name="log2_enrichment"),
            on=["cell_type", "cluster"],
            how="left",
            validate="one_to_one",
        )
    )
    source_data["global_mean_abundance"] = source_data["cell_type"].map(global_mean)
    source_data["n_spots"] = source_data["cluster"].map(cluster_spot_counts)
    source_data["n_samples"] = source_data["cluster"].map(cluster_sample_counts)
    source_data["plotted"] = source_data["cell_type"].isin(plotted_cell_types)
    source_data["abundance_key"] = abundance_key
    abundance_legend = size_legend_title or (
        "Mean q05 abundance\n(cells per spot)"
        if abundance_key == "q05_cell_abundance_w_sf"
        else "Mean posterior abundance\n(cells per spot)"
    )

    save_path = str(save_path)
    root, extension = os.path.splitext(save_path)
    if extension.lower() != ".pdf":
        save_path = f"{root if extension else save_path}.pdf"
        root = os.path.splitext(save_path)[0]
    if source_data_path is None:
        source_data_path = f"{root}_source.tsv"
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(source_data_path)), exist_ok=True)
    source_data.to_csv(source_data_path, sep="\t", index=False)

    n_rows, n_cols = plot_abundance.shape
    fig_width = min(7.2, max(4.6, 2.5 + 0.32 * n_cols))
    fig_height = min(8.5, max(3.2, 1.45 + 0.22 * n_rows))
    abundance_values = plot_abundance.to_numpy(dtype=float)
    positive_values = abundance_values[abundance_values > 0]
    size_cap = float(np.quantile(positive_values, 0.99)) if positive_values.size else 1.0
    size_cap = max(size_cap, np.finfo(float).eps)
    marker_sizes = 72.0 * np.clip(abundance_values / size_cap, 0, 1)

    xs, ys = np.meshgrid(np.arange(n_cols), np.arange(n_rows))
    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
        norm = TwoSlopeNorm(vmin=-enrichment_clip, vcenter=0, vmax=enrichment_clip)
        scatter = ax.scatter(
            xs.ravel(),
            ys.ravel(),
            c=np.clip(plot_enrichment.to_numpy(dtype=float), -enrichment_clip, enrichment_clip).ravel(),
            s=marker_sizes.ravel(),
            cmap="RdBu_r",
            norm=norm,
            edgecolors="#4D4D4D",
            linewidths=0.25,
            rasterized=n_rows * n_cols > 1500,
        )
        ax.set_xticks(np.arange(n_cols))
        ax.set_xticklabels(
            cluster_order,
            rotation=90 if n_cols > 15 else 45,
            ha="center" if n_cols > 15 else "right",
            fontsize=5.5 if n_cols > 15 else 6,
        )
        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(cell_type_order, fontsize=5.8)
        ax.invert_yaxis()
        ax.set_xlim(-0.55, n_cols - 0.45)
        ax.set_ylim(n_rows - 0.45, -0.55)
        ax.set_xlabel(x_axis_label, fontsize=7)
        ax.set_ylabel("Reference cell type", fontsize=7)
        title = plot_title or "Relative cell abundance across spatial regions"
        if omitted_cell_types:
            title += f" (top {len(plotted_cell_types)} of {len(ranked_cell_types)})"
        ax.set_title(title, loc="left", fontsize=8, fontweight="semibold", pad=5)
        ax.tick_params(axis="both", length=0, pad=2)
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_visible(False)

        cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.025, aspect=24)
        cbar.set_label("log2 enrichment\nvs tissue mean", fontsize=6.2, labelpad=4)
        cbar.ax.tick_params(labelsize=5.5, length=2, width=0.4)
        cbar.outline.set_linewidth(0.4)

        legend_values = np.unique(np.linspace(size_cap / 3, size_cap, 3))
        legend_handles = [
            ax.scatter(
                [], [],
                s=72.0 * min(value / size_cap, 1),
                facecolor="#BDBDBD",
                edgecolor="#4D4D4D",
                linewidth=0.25,
            )
            for value in legend_values
        ]
        legend_labels = [f"{value:.2g}" for value in legend_values]
        legend_labels[-1] = f"≥{legend_labels[-1]}"
        ax.legend(
            legend_handles,
            legend_labels,
            title=abundance_legend,
            loc="lower left",
            bbox_to_anchor=(1.01, 0.02),
            frameon=False,
            fontsize=5.5,
            title_fontsize=5.8,
            labelspacing=0.7,
            handletextpad=0.6,
            borderaxespad=0,
        )
        fig.savefig(
            save_path,
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.04,
            facecolor="white",
            metadata={"Title": plot_title or "Relative cell abundance across spatial regions"},
        )
        plt.close(fig)

    logger.info("%s regional abundance dotplot saved to %s", method_label, save_path)
    logger.info("%s dotplot source data saved to %s", method_label, source_data_path)
    return {
        "abundance_key": abundance_key,
        "abundance_matrix": abundance_df,
        "mean_abundance": mean_abundance,
        "log2_enrichment": log2_enrichment,
        "source_data": source_data,
        "plotted_cell_types": plotted_cell_types,
        "omitted_cell_types": omitted_cell_types,
        "cluster_spot_counts": cluster_spot_counts,
        "cluster_sample_counts": cluster_sample_counts,
        "pdf_path": save_path,
        "source_data_path": str(source_data_path),
    }
  
  
  
  
  
