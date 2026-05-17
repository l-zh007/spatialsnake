import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import os
import spatialdata as spd
import scanpy as sc
import scanpy.external as sce
import json
import gc
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import sparse
from scipy.stats import pearsonr, spearmanr

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
        print(f"图表已保存到: {save_path}")
    
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
        print(f"图表已保存到: {save_path}")
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
def plot_cell2location_dotplot(adata, unsupervised_cluster_col, save_path):
    def _to_dataframe(x, index=None, columns=None):
        if isinstance(x, pd.DataFrame):
            df = x.copy()
            if index is not None:
                df.index = index
            if columns is not None and len(columns) == df.shape[1]:
                df.columns = columns
            return df
        if sparse.issparse(x):
            x = x.toarray()
        x = np.asarray(x)
        if x.ndim != 2:
            raise ValueError(f"Expected a 2D abundance matrix, got shape={x.shape}")
        return pd.DataFrame(x, index=index, columns=columns)

    def _extract_cell2location_abundance(adata):
        factor_names = None
        mod_uns = adata.uns.get("mod", {})
        if isinstance(mod_uns, dict):
            factor_names = mod_uns.get("factor_names")

        candidate_keys = [
            "q05_cell_abundance_w_sf",
            "means_cell_abundance_w_sf",
            "q95_cell_abundance_w_sf",
            "stds_cell_abundance_w_sf",
        ]

        for key in candidate_keys:
            if key in adata.obsm:
                abundance = _to_dataframe(
                    adata.obsm[key],
                    index=adata.obs_names,
                    columns=factor_names,
                )
                if abundance.columns.isna().any() or abundance.columns.duplicated().any():
                    abundance.columns = [f"cell_type_{i}" for i in range(abundance.shape[1])]
                return abundance, key

        if factor_names is not None:
            prefix_candidates = [
                "q05_cell_abundance_w_sf_",
                "means_cell_abundance_w_sf_",
                "q95_cell_abundance_w_sf_",
            ]
            for prefix in prefix_candidates:
                cols = [f"{prefix}{ct}" for ct in factor_names if f"{prefix}{ct}" in adata.obs.columns]
                if cols:
                    abundance = adata.obs[cols].copy()
                    abundance.columns = [c.replace(prefix, "", 1) for c in cols]
                    abundance.index = adata.obs_names
                    return abundance, prefix.rstrip("_")

        raise KeyError(
            "No cell2location abundance matrix found. Expected one of the standard keys "
            "in adata.obsm (q05_cell_abundance_w_sf / means_cell_abundance_w_sf / q95_cell_abundance_w_sf) "
            "or column-wise abundance values in adata.obs."
        )

    def _correlate(abundance_df, cluster_series, method):
        cluster_dummies = pd.get_dummies(cluster_series.astype(str), prefix=None)
        corr = pd.DataFrame(index=abundance_df.columns, columns=cluster_dummies.columns, dtype=float)

        for cell_type in abundance_df.columns:
            x = abundance_df[cell_type].astype(float).to_numpy()
            for cluster in cluster_dummies.columns:
                y = cluster_dummies[cluster].astype(float).to_numpy()

                valid = np.isfinite(x) & np.isfinite(y)
                x_valid = x[valid]
                y_valid = y[valid]

                if len(x_valid) < 3 or np.std(x_valid) == 0 or np.std(y_valid) == 0:
                    value = np.nan
                elif method == "pearson":
                    value = pearsonr(x_valid, y_valid)[0]
                elif method == "spearman":
                    value = spearmanr(x_valid, y_valid)[0]
                else:
                    raise ValueError("method must be 'pearson' or 'spearman'")

                corr.loc[cell_type, cluster] = value

        return corr

    if unsupervised_cluster_col not in adata.obs.columns:
        raise KeyError(
            f"'{unsupervised_cluster_col}' not found in adata.obs. "
            f"Available columns include: {list(adata.obs.columns[:10])} ..."
        )

    abundance_df, abundance_key = _extract_cell2location_abundance(adata)
    cluster_series = adata.obs[unsupervised_cluster_col].copy()

    valid_mask = cluster_series.notna()
    abundance_df = abundance_df.loc[valid_mask]
    cluster_series = cluster_series.loc[valid_mask]

    abundance_df = abundance_df.loc[:, abundance_df.var(axis=0) > 0]
    if abundance_df.shape[1] == 0:
        raise ValueError("All abundance columns have zero variance after filtering.")

    methods = ["pearson", "spearman"]
    corr_results = {m: _correlate(abundance_df, cluster_series, m) for m in methods}

    plot_method = "pearson"
    corr_df = corr_results[plot_method].copy()

    ordered_celltypes = (
        corr_df.abs()
        .max(axis=1)
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    ordered_clusters = (
        corr_df.abs()
        .max(axis=0)
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    corr_df = corr_df.loc[ordered_celltypes, ordered_clusters]

    n_rows = len(ordered_celltypes)
    n_cols = len(ordered_clusters)

    fig_width = max(7.5, 0.95 * n_cols + 2.4)
    fig_height = max(7.5, 0.42 * n_rows + 1.8)

    fig, ax = plt.subplots(
        1,
        1,
        figsize=(fig_width, fig_height),
        dpi=320,
    )

    cmap = plt.cm.PuOr
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

    xs, ys, colors, sizes, labels = [], [], [], [], []
    for yi, cell_type in enumerate(corr_df.index):
        for xi, cluster in enumerate(corr_df.columns):
            value = corr_df.loc[cell_type, cluster]
            if pd.isna(value):
                continue
            xs.append(xi)
            ys.append(yi)
            colors.append(value)
            sizes.append(90 + 900 * abs(value))
            labels.append(f"{value:.2f}")

    scatter = ax.scatter(
        xs,
        ys,
        c=colors,
        s=sizes,
        cmap=cmap,
        norm=norm,
        edgecolors="white",
        linewidths=0.9,
        alpha=0.98,
    )

    for x, y, label, value in zip(xs, ys, labels, colors):
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=7.5,
            color="white" if abs(value) > 0.42 else "black",
            fontweight="bold",
        )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(corr_df.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(corr_df.index, fontsize=9)
    ax.invert_yaxis()

    ax.set_xlim(-0.55, n_cols - 0.45)
    ax.set_ylim(n_rows - 0.45, -0.55)

    ax.set_xlabel(unsupervised_cluster_col, fontsize=11, labelpad=8)
    ax.set_ylabel("Cell2location cell types", fontsize=11, labelpad=8)
    ax.set_title("Pearson correlation", fontsize=13, fontweight="bold", pad=8)

    ax.set_facecolor("#f7f7f7")
    ax.set_axisbelow(True)
    ax.grid(color="#d9d9d9", linestyle="--", linewidth=0.45, alpha=0.55)

    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(scatter, ax=ax, fraction=0.028, pad=0.015, aspect=28)
    cbar.set_label("r", fontsize=9, labelpad=6)
    cbar.ax.tick_params(labelsize=8, length=2)
    cbar.outline.set_linewidth(0.4)

    fig.suptitle(
        f"Cell2location abundance vs {unsupervised_cluster_col}",
        fontsize=15,
        fontweight="bold",
        y=0.995,
    )

    fig.subplots_adjust(
        left=0.23,
        right=0.92,
        bottom=0.14,
        top=0.95,
    )

    save_path = str(save_path)
    if not save_path.lower().endswith(".png"):
        save_path = f"{save_path}.png"

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    plt.savefig(save_path, dpi=320, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)

    return {
        "abundance_key": abundance_key,
        "abundance_matrix": abundance_df,
        "cluster_series": cluster_series,
        "correlation": corr_results,
        "plotted_method": plot_method,
    }
  
  
  
  
  
