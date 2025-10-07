import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

from scipy import sparse
import pytacs as tax
import scanpy as sc
import spatialdata as sd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
####'../xenium_multiple/breast_cancer.zarr'
zarr_save_path = './breast2.zarr'
sc_data_path = 'GSM7782698_count_raw_feature_bc_matrix.h5'
cell_type_path = 'GSE243275_Barcode_Cell_Type_Matrices.xlsx'


sp_adata = sd.read_zarr(zarr_save_path).tables["table"]
bcs = pd.read_excel(cell_type_path)
sc_adata = sc.read_10x_h5(sc_data_path)

print("空间数据形状:", sp_adata.var)
print("单细胞数据形状:", sc_adata)
print("注释数据形状:", bcs.shape)

sp_adata.var["gene_name"] = sp_adata.var.index
sc_adata.var["gene_name"] = sc_adata.var.index
sc_adata.var.set_index("gene_ids", inplace=True)
sc_adata.var["gene_ids"] = sc_adata.var.index

print("单细胞数据形状:", sc_adata.var)







sc_gene_ids = sc_adata.var.index.tolist()
sc_gene_ids_valid = [g for g in sc_gene_ids if not g.startswith("DEPRECATED_")]
print(f"sc_adata 总gene_ids数：{len(sc_gene_ids)}")
print(f"sc_adata 有效gene_ids数（过滤后）：{len(sc_gene_ids_valid)}")

sp_gene_ids = sp_adata.var['gene_ids'].tolist()
sp_gene_ids_valid = [g for g in sp_gene_ids if pd.notna(g) and str(g).strip() != ""]
print(f"sp_adata 总gene_ids数：{len(sp_gene_ids)}")
print(f"sp_adata 有效gene_ids数（过滤后）：{len(sp_gene_ids_valid)}")


sc_gene_set = set(sc_gene_ids_valid)
sp_gene_set = set(sp_gene_ids_valid)
common_gene_ids = sc_gene_set & sp_gene_set
common_gene_count = len(common_gene_ids)
common_gene_ids=list(common_gene_ids)
print(sp_adata.var['gene_ids'])



sc_adata = sc_adata[bcs['Barcode'], common_gene_ids]
sc_adata.obs = bcs
sc_adata.obs['celltype'] = sc_adata.obs['Annotation']

sp_adata = sp_adata[:, sp_adata.var['gene_ids'].isin(common_gene_ids)].copy()

sp_adata.var = sc_adata.var



# if 'gene_ids' not in sc_adata.var.columns:
#     sc_adata.var['gene_ids'] = sc_adata.var.index
# 
# if 'gene_ids' not in sp_adata.var.columns:
#     sp_adata.var['gene_ids'] = sp_adata.var.index

# shared_genes = np.intersect1d(sc_adata.var['gene_ids'], sp_adata.var['gene_ids'])


print(sc_adata)
print(sp_adata)
print(sc_adata.var)
print(sp_adata.var)
print(sc_adata.obs)
print(sp_adata.obs)



print(sc_adata.X)
print(sp_adata.X)



data_prep = tax.AnnDataPreparer(sc_adata, sp_adata,sn_colname_celltype="Annotation")
clf = tax.SVM()
clf.fit(data_prep.sn_adata)
agg_res = tax.rw_aggregate(
    st_anndata=data_prep.sp_adata,
    classifier=clf,
    max_iter=20,
    steps_per_iter=3,
    nbhd_radius=2.4,
    max_propagation_radius=10.,
    mode_metric='inv_dist',
    mode_embedding='pc',
    mode_aggregation='unweighted',
    n_pcs=50,
    verbose=False
)
ct_full = tax.extract_celltypes_full(agg_res)
print(agg_res.dataframe)



def plot_full_spatial_distribution(ct_full, spatial_coords):
    cell_types = ct_full
    unique_types = np.unique(cell_types)
    tab20 = plt.cm.get_cmap('tab20', 20)
    colors = tab20(np.arange(20))
    color_map = {}
    for i, cell_type in enumerate(unique_types):
        if cell_type == 'Undefined':
            color_map[cell_type] = [0.9, 0.9, 0.9, 1]
        else:
            color_map[cell_type] = colors[i % 20]
    cell_colors = np.array([color_map[ct] for ct in cell_types])
    plt.figure(figsize=(12, 8), dpi=300)
    plt.scatter(
        x=spatial_coords[:, 0],
        y=spatial_coords[:, 1],
        s=1,
        c=cell_colors,
        edgecolors='none'
    )
    legend_elements = []
    for cell_type in unique_types:
        if cell_type != 'Undefined':
            legend_elements.append(Line2D(
                [0], [0],
                marker='o',
                color='w',
                markerfacecolor=color_map[cell_type],
                markersize=10,
                label=cell_type
            ))
    legend_elements.append(Line2D(
        [0], [0],
        marker='o',
        color='w',
        markerfacecolor=[0.9, 0.9, 0.9],
        markersize=10,
        label='Undefined'
    ))
    plt.legend(
        handles=legend_elements,
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        ncol=1,
        frameon=False,
        title='Cell Types'
    )
    plt.title('Spatial Distribution of Cell Types', fontsize=24)
    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()
    plt.show()
    plt.savefig(
        "./TCA.png",
        dpi=300,
        bbox_inches='tight')
spatial_coords = data_prep.sp_adata.obsm['spatial']
plot_full_spatial_distribution(ct_full, spatial_coords)






