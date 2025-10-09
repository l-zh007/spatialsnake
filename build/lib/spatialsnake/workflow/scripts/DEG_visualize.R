library(tidyverse)
library(pheatmap)
library(RColorBrewer)
filename<-".//visium_multiple//DEG_results//aggregated_diffexp.csv"
aggregated_diffexp<-read_csv(filename)
colnames(aggregated_diffexp)[1] <-"gene_name"

aggregated_diffexp$significant <- ifelse(
  aggregated_diffexp$padj < 0.5 & abs(aggregated_diffexp$log2FoldChange) > 1,
  ifelse(aggregated_diffexp$log2FoldChange > 1, "Up", "Down"),
  "Not significant")

re_df<- filter(aggregated_diffexp,abs(log2FoldChange)>1 & pvalue < 0.05 & baseMean!=0)
re_df<-as.data.frame(re_df)



ggplot(aggregated_diffexp, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point() +
  # 自定义颜色映射，up 为红色，down 为蓝色
  scale_color_manual(values = c("Up" = "red", "Down" = "blue","Not significant"= "grey","NA"="grey")) +
  labs(
    title = "cancer and normal",
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "type"
  ) +
  theme_minimal()


# 绘制 MA 图
ggplot(re_df, aes(x = log2(re_df$mean_all), y = log2FoldChange, color = direction)) +
  geom_point() +
  # 自定义颜色：红色（上调）、蓝色（下调）、灰色（不显著）
  scale_color_manual(values = c("up" = "red", "down" = "blue", "not significant" = "gray")) +
  labs(
    title = "control_vs_AGES_LPS",  # 标题，按实际对比组修改
    x = "Log2 mean readcount",       # 横坐标标签
    y = "Log2 fold change",          # 纵坐标标签
    color = "Type"                   # 图例标题
  ) +
  theme_minimal()






library(pheatmap)
library(RColorBrewer)

res<-results(dds)
sig_genes <- rownames(res[res$padj < 0.5 & abs(res$log2FoldChange) > 1, ])


rld <- rlog(dds, blind = FALSE)  
expr_matrix <- assay(rld)[sig_genes, ] 



sample_names <- colnames(expr_matrix)  # 如 c("ABA_1", "ABA_2", "EtOH_1", "EtOH_2")
sample_annot <- data.frame(
  Group = factor(c("ABA", "ABA", "EtOH", "EtOH")), 
  row.names = sample_names)


color_pal <- colorRampPalette(c("blue", "white", "red"))(100)  

anno_colors <- list(
  Group = c(ABA = "#00BFC4", EtOH = "#F8766D")
)

# 绘制热图
pheatmap(
  expr_matrix,                # 标准化表达矩阵
  annotation_col = sample_annot,  # 样本分组注释
  annotation_colors = anno_colors, # 分组颜色映射
  annotation_names_col = FALSE,    # 隐藏注释列名
  show_rownames = FALSE,           # 隐藏基因名（基因多时间更清晰）
  show_colnames = TRUE,            # 显示样本名
  cluster_rows = TRUE,             # 基因聚类（揭示表达模式）
  cluster_cols = TRUE,             # 样本聚类（揭示重复一致性）
  color = color_pal,               # 颜色梯度
  main = "差异基因聚类热图",       # 标题
  fontsize = 10,                   # 字体大小
  border_color = NA,               # 去除单元格边框
  treeheight_row = 20,             # 基因聚类树高度
  treeheight_col = 20              # 样本聚类树高度
)





















