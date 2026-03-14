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

ensure_output_file <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    write.csv(data.frame(), path, row.names = FALSE)
  }
}

ensure_output_files <- function(paths) {
  for (p in paths) {
    ensure_output_file(p)
  }
}

ensure_output_file(opt$output_path)

get_pval_col <- function(df) {
  if ("padj" %in% colnames(df)) {
    return("padj")
  }
  if ("FDR" %in% colnames(df)) {
    return("FDR")
  }
  "pvalue"
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

prepare_significance <- function(df, pcol, p_cut, fc_cut) {
  df$pval_use <- df[[pcol]]
  df$significant <- ifelse(
    !is.na(df$pval_use) & df$pval_use < p_cut & abs(df$log2FoldChange) > fc_cut,
    ifelse(df$log2FoldChange > 0, "Up", "Down"),
    "Not significant"
  )
  df
}

if (opt$algorithm == "edgeR"){
  counts <- read.csv(opt$input_dir, row.names = 1, header = TRUE)
  counts<-column_to_rownames(counts,var='...1')
  group_path <- file.path(parent_dir, "group.csv")
  fallback_group_path <- file.path(parent_dir, "edgeR_counts.csv")
  if (file.exists(group_path)) {
    group <- read.csv(group_path, header = TRUE)
  } else {
    group <- read.csv(fallback_group_path, header = TRUE)
  }
  if ("condition" %in% colnames(group)) {
    group <- group$condition
  } else if ("group" %in% colnames(group)) {
    group <- group$group
  } else {
    group <- group[, 1]
  }
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
}

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


run_enrich <- function(gene_name, sapcies, path) {
  go_path <- file.path(path, "GO_data.csv")
  kegg_path <- file.path(path, "kegg_data.csv")
  ensure_output_files(c(go_path, kegg_path))
  if (length(gene_name) == 0) {
    return(NULL)
  }
  tryCatch({
    Gene_ID <- bitr(gene_name, fromType="SYMBOL",
                             toType="ENTREZID",
                             OrgDb="org.Hs.eg.db")
    if (is.null(Gene_ID) || nrow(Gene_ID) == 0) {
      return(NULL)
    }
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
    write.csv(GOenrich@result, go_path)
    write.csv(kegg_enrich@result, kegg_path)
    enrich_bubble_by_ontology(GOenrich, path)
    enrich_kegg_fun(kegg_enrich,path)
  }, error = function(e) {
    ensure_output_files(c(go_path, kegg_path))
  })
}

if (!"gene" %in% colnames(aggregated_diffexp)) {
  aggregated_diffexp$gene <- rownames(aggregated_diffexp)
}

contrast_list <- if ("contrast" %in% colnames(aggregated_diffexp)) unique(aggregated_diffexp$contrast) else c("comparison")
p_cut_strict <- 0.05
fc_cut_strict <- 1
p_cut_loose <- 0.1
fc_cut_loose <- 0.5
pval_col <- get_pval_col(aggregated_diffexp)

summary_rows <- list()

for (contrast_name in contrast_list) {
  if ("contrast" %in% colnames(aggregated_diffexp)) {
    df_ct <- aggregated_diffexp %>% filter(contrast == contrast_name)
  } else {
    df_ct <- aggregated_diffexp
  }
  if (nrow(df_ct) == 0) {
    next
  }
  df_ct <- prepare_significance(df_ct, pval_col, p_cut_strict, fc_cut_strict)
  df_ct_loose <- prepare_significance(df_ct, pval_col, p_cut_loose, fc_cut_loose)
  ct_dir <- if (contrast_name == contrast_list[1]) parent_dir else file.path(parent_dir, safe_name(contrast_name))
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df_ct, file.path(ct_dir, "diff_all.csv"))

  p9<-ggplot(df_ct, aes(x = log2FoldChange, y = -log10(pval_use), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue","Not significant"= "grey","NA"="grey")) +
    labs(
      title = contrast_name,
      x = "log2FoldChange",
      y = paste0("-log10(", pval_col, ")"),
      color = "type"
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(ct_dir, "vocanal.pdf"),
    plot = p9,
    width = 12,   
    height = 8,   
    units = "in")

  positive_strict <- filter(df_ct, significant == "Up")
  negative_strict <- filter(df_ct, significant == "Down")
  positive_loose <- filter(df_ct_loose, significant == "Up")
  negative_loose <- filter(df_ct_loose, significant == "Down")

  if (nrow(positive_loose) + nrow(negative_loose) == 0) {
    df_ct_loose <- prepare_significance(df_ct, pval_col, 0.2, 0.25)
    positive_loose <- filter(df_ct_loose, significant == "Up")
    negative_loose <- filter(df_ct_loose, significant == "Down")
  }

  positive_file_path <- file.path(ct_dir, "positive")
  negative_file_path <- file.path(ct_dir, "negative")
  dir.create(positive_file_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(negative_file_path, recursive = TRUE, showWarnings = FALSE)

  write.csv(positive_strict, file.path(positive_file_path, "diff_strict.csv"))
  write.csv(negative_strict, file.path(negative_file_path, "diff_strict.csv"))
  write.csv(positive_loose, file.path(positive_file_path, "diff_loose.csv"))
  write.csv(negative_loose, file.path(negative_file_path, "diff_loose.csv"))

  pos_genes <- rownames(positive_strict)
  neg_genes <- rownames(negative_strict)
  if (length(pos_genes) == 0) {
    pos_genes <- rownames(positive_loose)
  }
  if (length(neg_genes) == 0) {
    neg_genes <- rownames(negative_loose)
  }

  run_enrich(pos_genes, opt$spacies, positive_file_path)
  run_enrich(neg_genes, opt$spacies, negative_file_path)

  summary_rows[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    up_strict = nrow(positive_strict),
    down_strict = nrow(negative_strict),
    up_loose = nrow(positive_loose),
    down_loose = nrow(negative_loose)
  )
}

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(parent_dir, "contrast_summary.csv"), row.names = FALSE)

if ("contrast" %in% colnames(aggregated_diffexp)) {
  heat_df <- aggregated_diffexp %>%
    select(gene, contrast, log2FoldChange) %>%
    distinct(gene, contrast, .keep_all = TRUE)
  heat_wide <- heat_df %>%
    tidyr::pivot_wider(names_from = contrast, values_from = log2FoldChange)
  if (nrow(heat_wide) > 1) {
    heat_mat <- as.matrix(heat_wide[, -1, drop = FALSE])
    rownames(heat_mat) <- heat_wide$gene
    heat_mat[is.na(heat_mat)] <- 0
    max_abs <- apply(abs(heat_mat), 1, max)
    top_n <- min(50, length(max_abs))
    top_genes <- names(sort(max_abs, decreasing = TRUE))[1:top_n]
    heat_mat <- heat_mat[top_genes, , drop = FALSE]
    if (nrow(heat_mat) > 1 && ncol(heat_mat) > 1 && length(unique(as.vector(heat_mat))) > 1) {
      color_pal <- colorRampPalette(c("blue", "white", "red"))(100)
      pdf(file.path(parent_dir, "contrast_log2fc_heatmap.pdf"), width = 12, height = 8)
      pheatmap(
        heat_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = color_pal,
        border_color = NA,
        show_rownames = FALSE,
        show_colnames = TRUE,
        fontsize = 9
      )
      dev.off()
    }
  }
}





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





