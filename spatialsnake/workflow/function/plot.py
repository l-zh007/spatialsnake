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
        dendrogram=True,
        standard_scale='var',  # 按TF标准化
        swap_axes=True,  # 交换轴，使TF在y轴
        show_gene_labels=True)
    plt.savefig(outputs, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()  
  
  
  
  
  
  
  
  
  
  
  
  
