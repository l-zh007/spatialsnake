devtools::install_github('zktuong/ktplots', dependencies = TRUE)
library(tidyverse)
library(ktplots)
library(SingleCellExperiment)
library(SingleCellExperiment)
library(reticulate)
library(SingleCellExperiment)
library(reticulate)
library(ktplots)
library(zellkonverter)
h5ad_path <- "./visium_multiple/cellphone/Normal_P3.h5ad"

sce <- readH5AD(h5ad_path)

use_condaenv("snakemake_env", required = TRUE)
py_config()
ad=import('anndata')

adata = ad$read_h5ad(h5ad_path)
counts <- Matrix::t(adata$X)
row.names(counts) <- row.names(adata$var)
colnames(counts) <- row.names(adata$obs)
sce <- SingleCellExperiment(list(counts = counts), colData = adata$obs, rowData = adata$var)
colnames(colData(sce))[12] <- "celltype"




means <- read.delim('./visium_multiple/cellphone/cellphonedb_output_Normal/statistical_analysis_means_Normal.txt', check.names = FALSE)
pvalues <- read.delim('./visium_multiple/cellphone/cellphonedb_output_Normal/statistical_analysis_pvalues_Normal.txt', check.names = FALSE)
deconvoluted <- read.delim('./visium_multiple/cellphone/cellphonedb_output_Normal/statistical_analysis_deconvoluted_Normal.txt', check.names = FALSE)
interaction_grouping <- read.delim('interactions_groups.txt')


plot_cpdb(
  cell_type1 = 'Fibroblast', 
  cell_type2 = 'Enterocyte', 
  scdata = sce,
  celltype_key = 'celltype',
  means = means, 
  pvals = pvals)


plot_cpdb_heatmap(
  pvals, 
  cell_types = NULL,
  scale = "none",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  border_color = "white",
  fontsize_row = 11,
  fontsize_col = 11,
  family = "Arial",
  main = "",
  treeheight_col = 0,
  treeheight_row = 0,
  low_col = "dodgerblue4",
  mid_col = "peachpuff",
  high_col = "deeppink4",
  alpha = 0.05,
  return_tables = FALSE,
  symmetrical = TRUE,   #设置为FALSE会有聚类的效果
)


plot_cpdb3(cell_type1 = 'Enterocyte', cell_type2 = 'Fibroblast',
           scdata = sce,
           celltype_key="celltype",                
           means = means,
           pvals = pvals,
           deconvoluted = deconvoluted, # new options from here on specific to plot_cpdb3
           keep_significant_only = TRUE,
           standard_scale = TRUE,
           remove_self = TRUE
)



plot_cpdb2(scdata = sce,
           cell_type1 = 'Enterocyte', # same usage style as plot_cpdb
           cell_type2 = 'Fibroblast',
           celltype_key = 'celltype',
           means = means,
           pvals = pvals,
           deconvoluted = deconvoluted, # new options from here on specific to plot_cpdb2
           gene_symbol_mapping = 'index', # column name in rowData holding the actual gene symbols if the row names is ENSG Ids. Might be a bit buggy
           desiredInteractions = list(c('Enterocyte', 'Fibroblast')),
           node_group_colors = c("Fibroblast" = "#86bc86", "Enterocyte" = "#ff7f0e"),
           keep_significant_only = TRUE,
           standard_scale = TRUE,
           remove_self = TRUE)

geneDotPlot(scdata = sce, # object
            genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"),
            celltype_key = 'celltype',
            standard_scale = TRUE) 





# p1 <- ccc_number_heatmap3('./cellphonedb_output_Normal/statistical_analysis_pvalues_.txt',title= 'Ctrl') #ggplot对象
# p2 <- ccc_number_heatmap3('./cellphonedb_output_Cancer/statistical_analysis_pvalues_.txt',title= 'PD') #ggplot对象
# write.csv(p1[["data"]],'1_Ctrl_interaction.csv')
# write.csv(p2[["data"]],'2_PD_interaction.csv')
# p1 + p2
# ggsave('3_Interaction.pdf',width = 15.8, height = 6.35)
# 
# 
# 
# 
# 
# 
# 
# 
# p3 <- ccc_bubble2(pfile="./cellphonedb_output_Ctrl/statistical_analysis_pvalues_.txt",
#                   mfile="./cellphonedb_output_Ctrl/statistical_analysis_means_.txt",
#                   top_n_pairs = 10,
#                   target.use = NULL,
#                   sort_by = "neg_log10",       # 新增：排序依据（可选"neg_log10"或"means_exp_log2"）
#                   source.use = c('Microglia'))
# p3        
# ggsave('4_Ctrl_Microglia_as_source_dotplot.pdf',width = 9.27, height = 5.77)
# # 组合使用其他参数
# p4 <- ccc_bubble2(pfile="./cellphonedb_output_PD/statistical_analysis_pvalues_.txt",
#                   mfile="./cellphonedb_output_PD/statistical_analysis_means_.txt",
#                   top_n_pairs = 10,
#                   target.use = NULL,
#                   sort_by = "neg_log10",       # 新增：排序依据（可选"neg_log10"或"means_exp_log2"）
#                   source.use = c('Microglia'))
# p4    
# ggsave('4_PD_Microglia_as_source.pdf',width = 9.27, height = 5.77)
# 
# 
# 
# 
# 
# 
# 
# ccc_compare_update(group1.name = "Ctrl",group2.name = "PD",
#                    group1.pfile = "./cellphonedb_output_Ctrl/statistical_analysis_pvalues_.txt",
#                    group1.mfile="./cellphonedb_output_Ctrl/statistical_analysis_means_.txt",
#                    group2.pfile="./cellphonedb_output_PD/statistical_analysis_pvalues_.txt",
#                    group2.mfile="./cellphonedb_output_PD/statistical_analysis_means_.txt",
#                    p.threshold = 0.01,thre=1,
#                    cell.pair=c(
#                      paste0("Microglia|",c("Astrocytes","DaNs")),
#                      paste0(c("Microglia",'OPCs'),"|Endothelial cells")
#                    ))
# ggsave('5_Ctrl_vs_PD.pdf',width = 8,height = 6.37)
# 
# 
# 
# 
# 
# 
# 
# ccc_compare2_update(group1.name = "Ctrl",group2.name = "PD",
#                     group1.pfile = "./cellphonedb_output_Ctrl/statistical_analysis_pvalues_.txt",
#                     group1.mfile="./cellphonedb_output_Ctrl/statistical_analysis_means_.txt",
#                     group2.pfile="./cellphonedb_output_PD/statistical_analysis_pvalues_.txt",
#                     group2.mfile="./cellphonedb_output_PD/statistical_analysis_means_.txt",
#                     p.threshold = 0.1,thre=0.5,
#                     cell.pair=c(
#                       paste0("Microglia|",c("Endothelial cells"))),
#                     plot.width=20,plot.height=15,filename = "6_Microglia_EC_"
# )
# ccc_line(table.path="/Volumes/Seagate Basic/精品代码/69_CellphoneDB实战/6_Microglia_EC_Ctrl2PD.xlsx",
#          ligand.cell="Microglia",receptor.cell="Endothelial cells",
#          group1.name = "Ctrl",group2.name = "PD",#这五个参数和上一步对应
#          ligand.color="#4dbbd6",receptor.color="#90d1c1",
#          pt.size=6,
#          line.thre1=0.5,line.thre2=6,#line.thre1和上一步的"thre"参数一致，line.thre2可以用来调整线的粗细，值越大，线越细
#          file.name="7_Microglia_EC_",plot.width=25,plot.height=20)















