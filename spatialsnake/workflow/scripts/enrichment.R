library(AnnotationDbi)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(optparse)
library(ggplot2)
option_list <- list(
  make_option(c("--species")),
  make_option(c("--input_dir")),
  make_option(c("--output_path")),
  make_option(c("--sample_id")),
  make_option(c("--type")),
  make_option(c("--input_zarr")),
  make_option(c("--input_spec"), default=""))
opt <- parse_args(OptionParser(option_list=option_list))
options(timeout = 6000000)
parent_dir <- dirname(opt$input_dir)
print(parent_dir)
file_path <- opt$output_path
output_path<-dirname(file_path)
species=opt$species

if (species == "human") {
  orgdb <- org.Hs.eg.db
  kegg_organism <- "hsa"
} else if (species == "mouse") {
  orgdb <- org.Mm.eg.db
  kegg_organism <- "mmu"
} else {
  stop("Unsupported species value: ", species, ". Expected 'human' or 'mouse'.")
}

filename<-opt$input_dir
marker_genes_pval <- read.csv(opt$input_dir, header = T)
marker_genes_pval <- marker_genes_pval[is.finite(marker_genes_pval$logfoldchanges) & marker_genes_pval$logfoldchanges > 0, ]
# Read only the tested gene index, without loading the expression matrix.
python_code <- paste(
  "import sys,zarr",
  "from spatialsnake.workflow.function.stereoseq_selection import resolve_stereoseq_table_key",
  "tables=zarr.open_group(sys.argv[1],mode='r')['tables']",
  "keys=list(tables.group_keys())",
  "key=resolve_stereoseq_table_key(sys.argv[3],keys) if sys.argv[2] in ['stereoseq','StereoSeq','Stereo-seq'] else ('table' if 'table' in keys else keys[0])",
  "var=tables[key]['var']",
  "print('\\n'.join(map(str,var[var.attrs['_index']][:])))",
  sep=";"
)
tested_genes <- system2("python", c("-c", shQuote(python_code), shQuote(opt$input_zarr), shQuote(opt$type), shQuote(opt$input_spec)), stdout=TRUE)
if (!is.null(attr(tested_genes, "status"))) stop("Could not read the tested gene universe")
universe <- bitr(unique(tested_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)$ENTREZID
# gene_list<-marker_genes_pval["gene_name"]
# gene_symbols <- gene_list$gene_name



group <- data.frame(gene=marker_genes_pval$names,
                    group=marker_genes_pval$group)

Gene_ID <- bitr(unique(marker_genes_pval$names), fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb=orgdb)

if (nrow(Gene_ID) == 0) {
  stop("No valid SYMBOL keys were found for species='", species,
       "'. Please check whether the marker gene symbols match the selected species.")
}

markers_data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

if (nrow(markers_data) == 0) {
  stop("No genes remained after merging ENTREZID conversion results with input groups.")
}


data_GO <- compareCluster(
  ENTREZID~group, 
  data=markers_data, 
  fun="enrichGO",
  universe=unique(universe), 
  OrgDb=orgdb,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_kegg <- compareCluster(
  ENTREZID~group, 
  data=markers_data, 
  fun="enrichKEGG",
  universe=unique(universe),
  organism = kegg_organism,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)


# Enrichment supports cluster annotation; it does not measure pathway activity.
# Keep complete result tables, including when no terms pass the cutoffs.
empty_results <- data.frame(Cluster = character(), ID = character(), Description = character(),
                            GeneRatio = character(), BgRatio = character(), pvalue = numeric(),
                            p.adjust = numeric(), qvalue = numeric(), geneID = character(), Count = integer())
go_results <- if (is.null(data_GO)) transform(empty_results, ONTOLOGY = character()) else data_GO@compareClusterResult
kegg_results <- if (is.null(data_kegg)) empty_results else data_kegg@compareClusterResult
write.csv(go_results, file.path(parent_dir, "GO_data.csv"))
write.csv(kegg_results, file.path(parent_dir, "kegg_data.csv"))

wrap_terms <- function(x) {
  vapply(x, function(term) paste(strwrap(term, width = 48), collapse = "\n"), character(1))
}

# One plotting theme and significance scale for all enrichment figures.
enrichment_theme <- theme_classic(base_size = 11, base_family = "sans") +
  theme(axis.text = element_text(colour = "#333333"),
        axis.line = element_line(linewidth = 0.35, colour = "#555555"),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major.x = element_line(colour = "#EAEAEA", linewidth = 0.25),
        legend.position = "right", plot.title = element_text(face = "bold"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

# The overview selects top five terms per cluster, then displays every
# significant cluster-term pair for that union (blank does not mean zero).
for (ont in c("BP", "CC", "MF")) {
  plot_data <- if (nrow(go_results)) {
    subset(go_results, ONTOLOGY == ont & is.finite(p.adjust) & p.adjust < 0.05)
  } else data.frame()
  if (nrow(plot_data)) {
    selected <- plot_data %>% arrange(p.adjust, ID) %>%
      group_by(Cluster) %>% slice_head(n = 5) %>% ungroup()
    term_ids <- unique(selected$ID[order(selected$p.adjust, selected$ID)])
    plot_data <- plot_data[plot_data$ID %in% term_ids, ]
    plot_data$ID <- factor(plot_data$ID, levels = rev(term_ids))
    plot_data$Cluster <- factor(as.character(plot_data$Cluster),
                                levels = unique(as.character(marker_genes_pval$group)))
    plot_data$ratio <- vapply(strsplit(as.character(plot_data$GeneRatio), "/", fixed = TRUE),
                              function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    labels <- setNames(wrap_terms(plot_data$Description), as.character(plot_data$ID))
    p <- ggplot(plot_data, aes(x = Cluster, y = ID)) +
      geom_point(aes(size = ratio, colour = p.adjust)) +
      scale_y_discrete(labels = labels) +
      scale_x_discrete(drop = FALSE) +
      scale_size_area(max_size = 5, name = "Gene ratio", limits = c(0, max(plot_data$ratio))) +
      scale_colour_gradient(low = "#08519C", high = "#C6DBEF", limits = c(0, 0.05),
                            name = "BH-adjusted P") +
      labs(x = "Cluster", y = NULL, title = paste("GO", ont, "· cluster markers")) +
      enrichment_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    height <- max(4, 1.6 + 0.22 * sum(lengths(strsplit(unname(labels[!duplicated(names(labels))]), "\n"))))
    width <- max(8, 6 + 0.35 * nlevels(plot_data$Cluster))
  } else {
    p <- ggplot() + annotate("text", x = 0, y = 0,
                             label = paste("No significant GO", ont, "terms")) + theme_void()
    width <- 8
    height <- 4
  }
  ggsave(file.path(parent_dir, paste0(ont, "_GO_cluster.png")), p,
         width = width, height = height, units = "in", dpi = 300, limitsize = FALSE)
}

# Per-cluster GO retains the three-ontology gene-count bars; KEGG uses a dot plot.
# Keep the existing five terms per ontology for GO and all returned KEGG terms.
for (analysis in c("GO", "KEGG")) {
  results <- if (analysis == "GO") go_results else kegg_results
  if (!nrow(results)) next
  count_limit <- max(results$Count)
  for (cluster in unique(results$Cluster)) {
    plot_data <- results[results$Cluster == cluster, ]
    plot_data <- plot_data[is.finite(plot_data$p.adjust) & plot_data$p.adjust < 0.05, ]
    if (analysis == "GO") {
      plot_data <- plot_data %>% arrange(p.adjust, ID) %>%
        group_by(ONTOLOGY) %>% slice_head(n = 5) %>% ungroup()
    }
    if (!nrow(plot_data)) next
    plot_data <- plot_data[order(plot_data$p.adjust, plot_data$ID), ]
    plot_data$ID <- factor(plot_data$ID, levels = rev(unique(plot_data$ID)))
    # Derive fold enrichment from the reported ratios for version compatibility.
    ratio <- function(x) vapply(strsplit(as.character(x), "/", fixed = TRUE),
                                function(y) as.numeric(y[1]) / as.numeric(y[2]), numeric(1))
    plot_data$fold <- ratio(plot_data$GeneRatio) / ratio(plot_data$BgRatio)
    labels <- setNames(wrap_terms(plot_data$Description), as.character(plot_data$ID))
    if (analysis == "GO") {
      p <- ggplot(plot_data, aes(x = Count, y = ID, fill = ONTOLOGY)) +
        geom_col(width = 0.8) +
        scale_y_discrete(labels = labels) +
        scale_fill_brewer(palette = "Set2", limits = c("BP", "CC", "MF"), drop = FALSE) +
        facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
        labs(x = "Gene Count", y = "Enriched Terms",
             title = paste0("Cluster ", cluster, " · GO enrichment")) +
        enrichment_theme + theme(legend.position = "none")
    } else {
      p <- ggplot(plot_data, aes(x = fold, y = ID)) +
        geom_point(aes(size = Count, colour = p.adjust)) +
        scale_y_discrete(labels = labels) +
        scale_size_area(max_size = 5, name = "Marker genes", limits = c(0, count_limit),
                        breaks = function(x) unique(pmax(1, round(pretty(x, n = 3))))) +
        scale_colour_gradient(low = "#08519C", high = "#C6DBEF", limits = c(0, 0.05),
                              name = "BH-adjusted P") +
        labs(x = "Fold enrichment", y = NULL,
             title = paste0("Cluster ", cluster, " · ", analysis, " enrichment")) + enrichment_theme
    }
    height <- max(4, 1.8 + 0.24 * sum(lengths(strsplit(unname(labels), "\n"))))
    dir.create(file.path(parent_dir, cluster), showWarnings = FALSE, recursive = TRUE)
    filename <- if (analysis == "GO") "GO_cluster.pdf" else "kegg_cluster.pdf"
    ggsave(file.path(parent_dir, cluster, filename), p, width = 9, height = height,
           units = "in", limitsize = FALSE)
  }
}














# plot_data <- data_kegg@compareClusterResult %>% filter(Cluster == 0)
# p3<-ggplot(plot_data, aes(x = -log10(pvalue))) +
#   geom_point(data = plot_data,
#              aes(x = -log10(pvalue),
#                  y = Description,
#                  size = Count,
#                  color = p.adjust)) +
#   theme(
#     axis.text.y = element_blank(),
#   ) +
#   labs(x = "-log10(pvalue)",
#        y = "")
# 
# flow_data <- plot_data %>%
#   group_by(subcategory, Description) %>%
#   summarise(n = n(), .groups = "drop") 
# 
# sankey_data <- flow_data %>%
#   make_long(subcategory,Description,value = n)
# 
# p4<-ggplot(sankey_data, 
#        aes(x = x,                # 层级位置
#            next_x = next_x,      # 下一层级位置
#            node = node,          # 当前节点名称
#            next_node = next_node,# 下一层节点名称
#            fill = node,          # 节点填充色（按节点名称）
#            value = value,        # 流量大小（条带宽度）
#            label = node)) +      # 节点标签文字
#   geom_sankey(
#     flow.alpha = 0.5,     # 流量条带的透明度
#     flow.fill = "grey",   # 流量条带的填充色
#     flow.color = "grey80",# 流量条带的描边色
#     node.fill = "lightblue", # 节点的填充色（可自定义）
#     smooth = 8,           # 条带的平滑度
#     width = 0.08          # 节点的宽度
#   ) +
#   geom_sankey_text(
#     size = 3.2,    # 标签文字大小
#     color = "black"# 标签文字颜色
#   ) +
#   theme_void() +        # 去除默认坐标轴等元素，使图更简洁
#   theme(legend.position = "none") +  # 隐藏图例（若不需要）
#   labs(title = "Sankey Diagram: Subcategory → Description")
# 
# p5 <- p4 + theme(plot.margin = unit(c(0,17,0,0),units="cm"))
# 
# ggdraw() + draw_plot(p5) + draw_plot(p3, scale = 0.7, x = 0.35, y=-0.22, width=0.7, height=1.37)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(data_kegg@compareClusterResult,aes(y=Description,x=Count,fill=pvalue))+
#   geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
#   scale_fill_gradient(low = "red",high ="blue" )+
#   labs(title = "KEGG Pathways Enrichment",
#        x = "Gene number",
#        y = "Pathway")+
#   theme(axis.title.x = element_text(face = "bold",size = 16),
#         axis.title.y = element_text(face = "bold",size = 16),
#         legend.title = element_text(face = "bold",size = 16))+
#   theme_bw()
# 
# 
# 
# 
# 
# 
# 
# 
# dotplot(data_GO, x="GeneRatio",color="qvalue",showCategory = 3,split = "ONTOLOGY",font.size = 10,label_format = 30)+ 
#   facet_grid(ONTOLOGY~., scale='free')
# 
# dev.off()
# 
# pdf(file="GO柱状图.pdf",width = 6,height = 9)
# 
# barplot(data_GO@compareClusterResult, drop = TRUE, showCategory =5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# dev.off()
# 
# 
# 
# plot_data <- data_GO %>%
#   as.data.frame() %>%
#   filter(p.adjust < 0.05) %>%
#   arrange(p.adjust) %>%
#   group_by(Cluster) %>%
#   arrange(p.adjust, .by_group = TRUE) %>%
#   mutate(rank = row_number()) %>%
#   arrange(Cluster, rank) %>%
#   slice_head(n = 5) %>%
#   ungroup() %>%  
#   mutate(Description = factor(Description,levels = unique(Description[order(p.adjust)])))
# 
# 
# ggplot(plot_data, aes(x = Cluster, y = Description)) +
#   geom_point(aes(size = Count, color = p.adjust)) +
#   scale_color_gradientn(colors = c(c("#FF9999", "#FF6666", "#FF3333", "#CC0000"))) +  # 红色更显著
#   scale_size(range = c(3, 5)) +  # 气泡大小范围
#   labs(
#     x = "Cluster",
#     y = "GO Term",
#     color = "Adjusted p-value",
#     size = "Gene Count",
#     title = "Enrichment Bubble Plot"
#   ) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# ggsave(
#   filename = paste0(output_path,'/',opt$cluster,'_go_plot.png'),  # 注意参数名对应 --spatial-plot
#   plot = spatial_plot,
#   width = 6, height = 6,         # 空间图需更大高度展示细节
#   dpi = 300
# )
# 
# 
# 
# 
# 
# 
# ggplot(plot_data, aes(x = reorder(Description, -p.adjust), y = Count, fill = p.adjust)) +
#   geom_bar(stat = "identity") +
#   scale_fill_gradientn(colors = c("#FF9999", "#FF6666", "#FF3333", "#CC0000")) +  # 颜色表示显著性
#   coord_flip() +  # 旋转坐标轴，方便显示长名称
#   labs(
#     x = "KEGG",
#     y = "Number of Genes",
#     fill = "Adjusted p-value",
#     title = "Top 10 Significant Enriched Terms"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5),  # 标题居中
#     axis.text.y = element_text(size = 10)
#   )
# 
# ggsave(
#   filename = paste0(output_path,'/',opt$cluster,'_kegg_plot.png'),  # 注意参数名对应 --spatial-plot
#   plot = spatial_plot,
#   width = 6, height = 6,         # 空间图需更大高度展示细节
#   dpi = 300
# )









# data_melt$cluster <- as.numeric(as.character(data_melt$cluster))  
# data_melt <- markers %>% 
#   filter(p_val < 0.05) %>%
#   group_by(cluster) %>% 
#   arrange(desc(avg_log2FC), .by_group = TRUE) %>%
#   slice_head(n = 3) %>%
#   mutate(rank = row_number()) %>%
#   ungroup()
# data_melt <- data_melt %>% 
#   arrange(cluster, rank) %>%
#   mutate(gene = factor(gene, levels = unique(gene)))
# 
# levels(plot_data$Cluster)
# plot_data$Cluster<-as.numeric(plot_data$Cluster)





# dotplot(plot_data, showCategory=5,font.size = 8)
# data_GO_sim_fil <- data_GO_sim@compareClusterResult

# library(forcats)
# data_GO_sim_fil$Description <- as.factor(data_GO_sim_fil$Description)
# data_GO_sim_fil$Description <- fct_inorder(data_GO_sim_fil$Description)
# 
# ggplot(data_GO_sim_fil, aes(Cluster, Description)) +
#   geom_point(aes(fill=p.adjust, size=Count), shape=21)+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
#         axis.text = element_text(color = 'black', size = 5))+
#   scale_fill_gradient(low="purple",high="yellow")+
#   labs(x=NULL,y=NULL)+
#   coord_flip()
















# file_path <- opt$output_path
# output_path<-dirname(file_path)
# 
# species="human"
# orgdb=ifelse(species == "human", "org.Hs.0eg.db", "org.Mm.eg.db")
# filename<-opt$input_dir
# re_df<-read_csv(filename)
# 
# gene_list<-re_df["gene_name"]
# 
# gene_symbols <- gene_list$gene_name
# 
# gene_id_convert <- bitr(
#   geneID = gene_symbols,
#   fromType = "SYMBOL",
#   toType = "ENTREZID",
#   OrgDb = org.Hs.eg.db
# )
# head(gene_id_convert)
# 
# gene_list_entrez <- gene_id_convert$ENTREZID
# gene_list_entrez <- gene_id_convert$SYMBOL
# GOenrich <- enrichGO(gene = gene_list_entrez,
#                      OrgDb = "org.Hs.eg.db",
#                      keyType = "ENTREZID",
#                      ont = opt$GO_ont,
#                      pAdjustMethod = "BH",
#                      qvalueCutoff = 0.05)
# 
# sorted_ego <- ego[order(ego$p.adjust), ]
# top10_ego <- head(sorted_ego, 10)
# go_ids <- top10_ego$ID
# 
# write.csv(GOenrich@result, 
#           file = paste0(output_path,"Go_enrichment_results.csv"), 
#           row.names = FALSE)
# 
# 
# 
# 
# gene_list_entrez <- gene_id_convert$ENTREZID
# kegg_enrich <- enrichKEGG(gene = gene_list_entrez,
#                  organism = 'hsa',
#                  pAdjustMethod = "BH",
#                  qvalueCutoff = 0.05)
# 
# write.csv(kegg_enrich@result, paste0(output_path,"kegg_enrichment_results.csv"), row.names = FALSE)
# 
# 
# 
# # 筛选显著条目（例如p.adjust < 0.05）并取前10个
# plot_data <- kegg_enrich %>%
#   filter(p.adjust < 0.05) %>%  # 筛选显著结果
#   arrange(p.adjust) %>%       # 按调整后p值排序
#   head(10) %>%                # 取前10个
#   mutate(Description = str_wrap(Description, width = 30))  # 换行显示长名称
# 
# # 绘制条形图
# ggplot(plot_data, aes(x = reorder(Description, -p.adjust), y = Count, fill = p.adjust)) +
#   geom_bar(stat = "identity") +
#   scale_fill_gradientn(colors = c("#FF9999", "#FF6666", "#FF3333", "#CC0000")) +  # 颜色表示显著性
#   coord_flip() +  # 旋转坐标轴，方便显示长名称
#   labs(
#     x = "KEGG",
#     y = "Number of Genes",
#     fill = "Adjusted p-value",
#     title = "Top 10 Significant Enriched Terms"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5),  # 标题居中
#     axis.text.y = element_text(size = 10)
#   )
# 
# ggsave(
#   filename = paste0(output_path,'/',opt$cluster,'_kegg_plot.png'),  # 注意参数名对应 --spatial-plot
#   plot = spatial_plot,
#   width = 6, height = 6,         # 空间图需更大高度展示细节
#   dpi = 300
# )
# 
# 
# 
# 
# 
# 
# 
# plot_data <- GOenrich %>%
#   filter(p.adjust < 0.05) %>%
#   head(15) %>%
#   mutate(
#     GeneRatio = Count / as.numeric(sub("/.*", "", BgRatio)),
#     Description = str_wrap(Description, width = 30)
#   )
# 
# 
# ggplot(plot_data, aes(x = Cluster, y = Description)) +
#   geom_point(aes(size = Count, color = p.adjust)) +
#   scale_color_gradientn(colors = c(c("#FF9999", "#FF6666", "#FF3333", "#CC0000"))) +  # 红色更显著
#   scale_size(range = c(3, 5)) +  # 气泡大小范围
#   labs(
#     x = "Cluster",
#     y = "GO Term",
#     color = "Adjusted p-value",
#     size = "Gene Count",
#     title = "Enrichment Bubble Plot"
#   ) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# ggsave(
#   filename = paste0(output_path,'/',opt$cluster,'_go_plot.png'),  # 注意参数名对应 --spatial-plot
#   plot = spatial_plot,
#   width = 6, height = 6,         # 空间图需更大高度展示细节
#   dpi = 300
# )







# gene_id_express <- re_df$log2FoldChange
# names(gene_id_express) <- gene_id_convert$ENTREZID
# gene_id_express <- gene_id_express[!is.na(names(gene_id_express))]
# gene_id_express <- sort(gene_id_express, decreasing = TRUE)
# 
# gse_GO <- gseGO(geneList = gene_id_express,
#                 OrgDb = 'org.Hs.eg.db',
#                 pvalueCutoff = 0.1,
#                 pAdjustMethod = 'BH',
#                 keyType = 'ENTREZID')
# write.csv(gse_GO@result, paste0(output_path,"kegg_enrichment_results.csv"), row.names = FALSE)
# gseaplot(
#   x = gse_GO,
#   geneSetID = 'GO:0044248',
#   title = "GSEA Enrichment Plot",
#   color.line = "blue",
#   color.vline = "#FF0000",
#   gene_geom = "bar"  # 新增参数，指定基因标记以条形显示
# )
# 
# ggsave(
#   filename = opt$`spatial_plot`,  # 注意参数名对应 --spatial-plot
#   plot = spatial_plot,
#   width = 6, height = 6,         # 空间图需更大高度展示细节
#   dpi = 300
# )






# gse_kegg<-gseKEGG(geneList = gene_id_express,
#         OrgDb = 'org.Hs.eg.db',
#         pvalueCutoff = 0.1,
#         pAdjustMethod = 'BH',
#         keyType = 'ENTREZID')


# gseaplot(
#   x = gse_kegg,
#   geneSetID = 'GO:0044248',
#   title = "GSEA Enrichment Plot",
#   color.line = "blue",
#   color.vline = "#FF0000",
#   gene_geom = "bar"  # 新增参数，指定基因标记以条形显示
# )
# 
# 
# 
# 
# ggsave(
#   filename = opt$`spatial_plot`,  # 注意参数名对应 --spatial-plot
#   plot = spatial_plot,
#   width = 6, height = 6,         # 空间图需更大高度展示细节
#   dpi = 300
# )













# # 获取通路信息
# pathway_info <- lapply(kk$ID, function(x) {
#   tryCatch({
#     kegg_get(x, organism = 'hsa')
#   }, error = function(e) {
#     NULL
#   })
# })
# names(pathway_info) <- kk$ID

# # 提取一级通路信息并添加到富集结果中
# kk$top_level_pathway <- sapply(pathway_info, function(x) {
#   if (!is.null(x)) {
#     return(x$category[1])
#   } else {
#     return(NA)
#   }
# })

# library(org.At.tair.db)
# library(AnnotationDbi)
# library(dplyr)
# 
# # 假设 allgene 是你的拟南芥基因列表（TAIR 格式）
# allgene <- c("AT1G01010", "AT1G01020", "AT1G01030")  # 示例基因
# 
# # 获取 GO 注释
# go_annot <- AnnotationDbi::select(
#   org.At.tair.db,
#   keys = allgene,
#   keytype = "TAIR",
#   columns = c("GO", "ONTOLOGY")
# ) %>%
#   rename(term = GO, GO_domain = ONTOLOGY) %>%
#   filter(!is.na(term), !is.na(GO_domain))  # 过滤缺失值
# 
# # 统计每个 GO 条目的基因数目
# go_summary <- go_annot %>%
#   group_by(GO_domain, term) %>%
#   summarise(gene_number = n()) %>%
#   ungroup()
# 
# # 查看结果
# head(go_summary)
# 
# 
# 
# 
# df_top20$Description <- factor(df_top20$Description, levels = df_top20$Description)
# ggplot(df_top20, aes(x = Description, y = Count, fill = ONTOLOGY)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.6) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
#     axis.text.y = element_text(size = 8),
#     axis.title = element_text(size = 10),
#     legend.text = element_text(size = 8),
#     legend.title = element_text(size = 10),
#     plot.title = element_text(size = 12),
#     plot.margin = margin(20, 20, 20, 20)
#   ) +
#   labs(
#     x = "",
#     y = "Gene Number",
#     title = "Enriched GO Terms (c_vs_t)"
#   )
