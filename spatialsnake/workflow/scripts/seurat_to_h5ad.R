library(schard)
# 单细胞转换
snhx = schard::h5ad2seurat("../st_project/lzh_project/Spatialsnake/data/Human_Lymph_Node.h5ad")

# 单样本空转转换
visx = schard::h5ad2seurat_spatial("../test/results/useful_results/cluster_Lesional_1.h5ad")

# 多样本空转转换
visl = schard::h5ad2seurat_spatial("../test/results/useful_results/concatenated_sdata.h5ad",simplify = FALSE)












