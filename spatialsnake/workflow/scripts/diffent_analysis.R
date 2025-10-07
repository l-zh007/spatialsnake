library(AnnotationDbi)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(optparse)
library(ggplot2)
library(ggsankey)
library(cowplot)
library(edgeR)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
# install.packages("cowplot")
option_list <- list(
  make_option(c("--spacies")),
  make_option(c("--input_dir")),
  make_option(c("--output_path")),
  make_option(c("--type")),
  make_option(c("--algorithm")))
opt <- parse_args(OptionParser(option_list=option_list))
parent_dir <- dirname(opt$input_dir)

if (opt$algorithm == "edgeR"){
  counts <- read.csv(opt$input_dir, row.names = 1, header = TRUE)
  counts<-column_to_rownames(counts,var='...1')
  group<- read.csv(paste0(parent_dir,"/edgeR_counts.csv"), header = TRUE,)
  print(counts)
  print(group)
  group<-group$condition
  dge <- DGEList(counts = counts, group = group)
  keep <- rowSums(cpm(dge) > 1) >= 1
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dge, method = 'TMM')  
  bcv <- 0.1
  et <- exactTest(dge, dispersion = bcv^2)
  results <- topTags(et, n = nrow(dge$counts))
  aggregated_diffexp<-results$table
  colnames(aggregated_diffexp)<-c("log2FoldChange",'logCPM','pvalue','FDR')
}
if (opt$algorithm == "DEseq2"){
  aggregated_diffexp <- read.csv(opt$input_dir, row.names = 1, header = TRUE)
  print(head(aggregated_diffexp))
}

aggregated_diffexp$significant <- ifelse(
  aggregated_diffexp$pvalue < 0.5 & abs(aggregated_diffexp$log2FoldChange) > 1,
  ifelse(aggregated_diffexp$log2FoldChange > 1, "Up", "Down"),
  "Not significant")

p9<-ggplot(aggregated_diffexp, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue","Not significant"= "grey","NA"="grey")) +
  labs(
    title = "cancer and normal",
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "type"
  ) +
  theme_minimal()
ggsave(
  filename = paste0(parent_dir,"/vocanal.pdf"),
  plot = p9,
  width = 12,   
  height = 8,   
  units = "in")

positive<-filter(aggregated_diffexp,significant=="Up")
negative<-filter(aggregated_diffexp,significant=="Down")

print(head(positive))
print(head(negative))

enrich_bubble_by_ontology <- function(data_GO, path) {
    ontology_list <- unique(data_GO$ONTOLOGY)
    plot_data <- data_GO %>%
      as.data.frame() %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      group_by(ONTOLOGY) %>%
      arrange(pvalue, .by_group = TRUE) %>%
      mutate(rank = row_number()) %>%
      arrange(ONTOLOGY, rank) %>%
      slice_head(n = 5) %>%
      ungroup() %>%  
      mutate(Description = factor(Description,levels = unique(Description[order(p.adjust)])))
    p1<-ggplot(plot_data, aes(x = Count, y = Description, fill = ONTOLOGY)) +
      geom_col(width = 0.8) +
      facet_grid(ONTOLOGY ~ ., scales = "free_y") +
      scale_fill_brewer(palette = "Set2") +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 10),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position = "none") +
      labs(x = "Gene Count", y = "Enriched Terms")
    ggsave(
      filename = paste0(path,"/GO_enrich.pdf"),
      plot = p1,
      width = 12,   
      height = 8,   
      units = "in")
  }

# enrich_bubble_by_ontology(GOenrich@result, parent_dir)
# 
# barplot(kegg_enrich,showCategory = 30,title = 'KEGG Pathway')
# dotplot(kegg_enrich)

enrich_kegg_fun <- function(data_kegg,path) {
    plot_data <- data_kegg@result
    plot_data <- plot_data %>%
      as.data.frame() %>%
      arrange(pvalue) %>%
      slice_head(n = 15) %>%
      mutate(Description = factor(Description,levels = unique(Description[order(p.adjust)])))
    
    
    p3<-ggplot(plot_data, aes(x = pvalue)) +
      geom_point(data = plot_data,
                 aes(x = pvalue,
                     y = subcategory,
                     size = Count,
                     color = p.adjust)) +
      theme(
        axis.text.y = element_blank(),
      ) +
      labs(x = "FoldEnrichment",
           y = "") +
      scale_color_gradientn(
        colors = c("#FFCCCC", "#FF6666", "#CC0000"),  
        name = "p.adjust"
      )
    
    flow_data <- plot_data %>%
      group_by(subcategory, Description) %>%
      summarise(n = n(), .groups = "drop") 
    
    sankey_data <- flow_data %>%
      make_long(subcategory,Description,value = n)
    
    all_nodes <- unique(c(sankey_data$node, sankey_data$next_node))
    node_colors <- hcl.colors(length(all_nodes), palette = "Set2")  
    names(node_colors) <- all_nodes
    
    p4<-ggplot(sankey_data, 
               aes(x = x,                # 层级位置
                   next_x = next_x,      # 下一层级位置
                   node = node,          # 当前节点名称
                   next_node = next_node,# 下一层节点名称
                   fill = node,          # 节点填充色（按节点名称）
                   value = value,        # 流量大小（条带宽度）
                   label = node)) +      # 节点标签文字
      geom_sankey(
        flow.alpha = 0.5,     # 流量条带的透明度
        flow.fill = "grey",   # 流量条带的填充色
        flow.color = "grey80",# 流量条带的描边色
        smooth = 8,           # 条带的平滑度
        width = 0.08
      ) +
      scale_fill_manual(values = node_colors) +
      geom_sankey_text(
        size = 3.2,    # 标签文字大小
        color = "black"# 标签文字颜色
      ) +
      theme_void() +        # 去除默认坐标轴等元素，使图更简洁
      theme(legend.position = "none") +  # 隐藏图例（若不需要）
      labs(title = "Sankey Diagram: Subcategory → Description")
    
    p5 <- p4 + theme(plot.margin = unit(c(0,19,0,0),units="cm"))
    p6<-ggdraw() + draw_plot(p5) + draw_plot(p3, scale = 0.7, x = 0.35, y=-0.22, width=0.7, height=1.37)
    ggsave(
      filename = paste0(path, "/kegg_cluster.pdf"),
      plot = p6,
      width = 12,   
      height = 8,
      units = "in")
}


enrich_funct<- function(gene_name,sapcies,path) {
  Gene_ID <- bitr(gene_name, fromType="SYMBOL",
                           toType="ENTREZID",
                           OrgDb="org.Hs.eg.db")
  GOenrich <- enrichGO(gene = Gene_ID$ENTREZID,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05)
  options(timeout = 180)
  kegg_enrich <- enrichKEGG(gene = Gene_ID$ENTREZID,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05)
  write.csv(GOenrich@result,paste0(path,'/GO_data.csv'))
  write.csv(kegg_enrich@result,paste0(path,'/kegg_data.csv'))
  enrich_bubble_by_ontology(GOenrich, path)
  enrich_kegg_fun(kegg_enrich,path)
}
print(positive$gene_name)

positive_file_path <- file.path(parent_dir,"positive")
negative_file_path <- file.path(parent_dir, "negative")
dir.create(positive_file_path, recursive = TRUE, showWarnings = FALSE)
dir.create(negative_file_path, recursive = TRUE, showWarnings = FALSE)
enrich_funct(rownames(positive),opt$spacies,positive_file_path)
enrich_funct(rownames(negative),opt$spacies,negative_file_path)





# res<-results(dds)
# sig_genes <- rownames(res[res$padj < 0.5 & abs(res$log2FoldChange) > 1, ])
# 
# rld <- rlog(dds, blind = FALSE)  
# expr_matrix <- assay(rld)[sig_genes, ] 
# 
# 
# 
# sample_names <- colnames(expr_matrix)  # 如 c("ABA_1", "ABA_2", "EtOH_1", "EtOH_2")
# sample_annot <- data.frame(
#   Group = factor(c("ABA", "ABA", "EtOH", "EtOH")), 
#   row.names = sample_names)
# 
# 
# color_pal <- colorRampPalette(c("blue", "white", "red"))(100)  
# 
# anno_colors <- list(
#   Group = c(ABA = "#00BFC4", EtOH = "#F8766D")
# )
# 
# # 绘制热图
# pheatmap(
#   expr_matrix,                # 标准化表达矩阵
#   annotation_col = sample_annot,
#   annotation_colors = anno_colors,
#   annotation_names_col = FALSE,
#   show_rownames = FALSE,
#   show_colnames = TRUE,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   color = color_pal,
#   main = "差异基因聚类热图",
#   fontsize = 10,
#   border_color = NA,
#   treeheight_row = 20,
#   treeheight_col = 20 
# )









