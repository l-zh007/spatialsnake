library(AnnotationDbi)
library(clusterProfiler)
library(dplyr)
library(edgeR)
library(ggplot2)
library(optparse)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(tidyr)

option_list <- list(
  make_option(c("--spacies")),
  make_option(c("--input_dir")),
  make_option(c("--output_path")),
  make_option(c("--type")),
  make_option(c("--algorithm")),
  make_option(c("--sample_list"), type = "character", default = ""),
  make_option(c("--cut_off_pvalue"), type = "double", default = 0.05),
  make_option(c("--cut_off_logFC"), type = "double", default = 1.5)
)
opt <- parse_args(OptionParser(option_list = option_list))
parent_dir <- dirname(opt$input_dir)
dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
annotation_key_cache <- new.env(parent = emptyenv())

print(opt$sample_list)

write_csv_file <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
}

save_png_plot <- function(plot_obj, file_path, width = 10, height = 7, dpi = 320) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

wrap_text <- function(x, width = 36) {
  vapply(
    x,
    function(value) paste(strwrap(as.character(value), width = width), collapse = "\n"),
    character(1)
  )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

read_groups_from_sample_list <- function(sample_list_path) {
  if (is.null(sample_list_path) || sample_list_path == "" || !file.exists(sample_list_path)) {
    return(character(0))
  }
  sample_df <- try(read.table(sample_list_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (inherits(sample_df, "try-error") || ncol(sample_df) < 3) {
    return(character(0))
  }
  unique(as.character(sample_df[[3]]))
}

resolve_group_names <- function(df, contrast_name, sample_list_path = "") {
  group1_name <- NA_character_
  group2_name <- NA_character_
  if ("group1" %in% colnames(df) && "group2" %in% colnames(df)) {
    group1_values <- unique(na.omit(as.character(df$group1)))
    group2_values <- unique(na.omit(as.character(df$group2)))
    if (length(group1_values) >= 1) {
      group1_name <- group1_values[1]
    }
    if (length(group2_values) >= 1) {
      group2_name <- group2_values[1]
    }
  }
  if (is.na(group1_name) || is.na(group2_name) || tolower(group1_name) == "group1" || tolower(group2_name) == "group2") {
    contrast_parts <- strsplit(as.character(contrast_name), "_vs_", fixed = TRUE)[[1]]
    if (length(contrast_parts) == 2) {
      group1_name <- contrast_parts[1]
      group2_name <- contrast_parts[2]
    }
  }
  if (is.na(group1_name) || is.na(group2_name) || tolower(group1_name) == "group1" || tolower(group2_name) == "group2") {
    sample_groups <- read_groups_from_sample_list(sample_list_path)
    if (length(sample_groups) >= 2) {
      group1_name <- sample_groups[1]
      group2_name <- sample_groups[2]
    }
  }
  if (is.na(group1_name) || group1_name == "") {
    group1_name <- "group1"
  }
  if (is.na(group2_name) || group2_name == "") {
    group2_name <- "group2"
  }
  list(group1 = group1_name, group2 = group2_name)
}

add_group_display_columns <- function(df, group1_name, group2_name) {
  df$display_group <- ifelse(
    df$significant == "Up",
    group1_name,
    ifelse(df$significant == "Down", group2_name, "None")
  )
  df$display_group <- factor(df$display_group, levels = c(group2_name, "None", group1_name))
  df
}

blank_plot <- function(title_text, subtitle_text = "") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = title_text, size = 6, fontface = "bold") +
    annotate("text", x = 0, y = -0.18, label = subtitle_text, size = 4.2, colour = "#4B5563") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void()
}

get_species_config <- function(species) {
  species_lower <- tolower(trimws(as.character(species)))
  if (species_lower %in% c("mouse", "mus", "mm", "mmu")) {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Mouse enrichment requires the org.Mm.eg.db package.")
    }
    return(list(
      orgdb = get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db")),
      kegg = "mmu"
    ))
  }
  list(orgdb = org.Hs.eg.db, kegg = "hsa")
}

get_pval_col <- function(df) {
  if ("padj" %in% colnames(df)) {
    return("padj")
  }
  if ("FDR" %in% colnames(df)) {
    return("FDR")
  }
  if ("pvalue" %in% colnames(df)) {
    return("pvalue")
  }
  stop("No p-value column found. Expected one of: padj, FDR, pvalue.")
}

get_abundance_col <- function(df) {
  if ("baseMean" %in% colnames(df)) {
    return("baseMean")
  }
  if ("logCPM" %in% colnames(df)) {
    return("logCPM")
  }
  NULL
}

ensure_gene_column <- function(df) {
  if (!"gene" %in% colnames(df)) {
    df$gene <- rownames(df)
  }
  df$gene <- as.character(df$gene)
  df
}

prepare_diff_table <- function(df, p_cut, fc_cut) {
  df <- ensure_gene_column(df)
  if (!"log2FoldChange" %in% colnames(df)) {
    stop("The differential result table must contain a log2FoldChange column.")
  }
  pval_col <- get_pval_col(df)
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  df$pval_use <- as.numeric(df[[pval_col]])
  df <- df %>%
    filter(!is.na(gene), gene != "", !is.na(log2FoldChange)) %>%
    distinct(gene, .keep_all = TRUE)
  df$plot_p <- ifelse(is.na(df$pval_use) | df$pval_use <= 0, 1e-300, df$pval_use)
  df$neg_log10_p <- -log10(df$plot_p)
  df$significant <- ifelse(
    !is.na(df$pval_use) & df$pval_use < p_cut & df$log2FoldChange >= fc_cut,
    "Up",
    ifelse(
      !is.na(df$pval_use) & df$pval_use < p_cut & df$log2FoldChange <= -fc_cut,
      "Down",
      "Not significant"
    )
  )
  df$significant <- factor(df$significant, levels = c("Up", "Down", "Not significant"))
  df
}

load_edger_results <- function(input_path) {
  counts <- read.csv(input_path, header = TRUE, check.names = FALSE)
  if ("...1" %in% colnames(counts)) {
    rownames(counts) <- counts$...1
    counts$...1 <- NULL
  } else {
    rownames(counts) <- counts[[1]]
    counts[[1]] <- NULL
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  group_path <- file.path(parent_dir, "group.csv")
  if (!file.exists(group_path)) {
    stop("edgeR requires group.csv next to the input count matrix.")
  }
  group_df <- read.csv(group_path, header = TRUE, check.names = FALSE)
  group <- if ("condition" %in% colnames(group_df)) {
    group_df$condition
  } else if ("group" %in% colnames(group_df)) {
    group_df$group
  } else {
    group_df[[1]]
  }
  if (length(group) != ncol(counts)) {
    stop("The number of groups does not match the number of samples in the edgeR count matrix.")
  }
  dge <- DGEList(counts = counts, group = group)
  keep <- rowSums(cpm(dge) > 1) >= 1
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  et <- exactTest(dge, dispersion = 0.1^2)
  result_df <- topTags(et, n = nrow(dge$counts))$table
  result_df$gene <- rownames(result_df)
  colnames(result_df) <- c("log2FoldChange", "logCPM", "pvalue", "FDR", "gene")
  result_df
}

load_diff_results <- function(input_path, algorithm) {
  algorithm_lower <- tolower(as.character(algorithm))
  if (algorithm_lower == "edger") {
    return(load_edger_results(input_path))
  }
  if (algorithm_lower == "deseq2") {
    return(read.csv(input_path, row.names = 1, header = TRUE, check.names = FALSE))
  }
  stop("Unsupported algorithm: ", algorithm)
}

make_volcano_plot <- function(df, contrast_name, p_cut, fc_cut, group1_name, group2_name) {
  volcano_df <- df
  volcano_p_col <- if ("pvalue" %in% colnames(volcano_df)) "pvalue" else "pval_use"
  volcano_df$volcano_pvalue <- as.numeric(volcano_df[[volcano_p_col]])
  volcano_df$volcano_pvalue[is.na(volcano_df$volcano_pvalue) | volcano_df$volcano_pvalue <= 0] <- 1e-300
  volcano_df$Sig <- ifelse(
    volcano_df$volcano_pvalue < p_cut & abs(volcano_df$log2FoldChange) >= fc_cut,
    ifelse(volcano_df$log2FoldChange > fc_cut, group1_name, group2_name),
    "None"
  )
  volcano_df$Sig <- factor(volcano_df$Sig, levels = c(group2_name, "None", group1_name))
  volcano_df$neg_log10_pvalue <- -log10(volcano_df$volcano_pvalue)
  if ("baseMean" %in% colnames(volcano_df)) {
    label_df <- volcano_df %>%
      filter(Sig %in% c(group1_name, group2_name)) %>%
      arrange(desc(baseMean), volcano_pvalue, desc(abs(log2FoldChange))) %>%
      slice_head(n = 12)
  } else {
    label_df <- volcano_df %>%
      filter(Sig %in% c(group1_name, group2_name)) %>%
      arrange(volcano_pvalue, desc(abs(log2FoldChange))) %>%
      slice_head(n = 12)
  }
  ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_pvalue, colour = Sig)) +
    geom_point(alpha = 0.4, size = 3.5) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(p_cut), lty = 4, col = "black", lwd = 0.8) +
    geom_text(
      data = label_df,
      aes(label = gene),
      size = 3.2,
      vjust = -0.5,
      show.legend = FALSE,
      check_overlap = TRUE
    ) +
    scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) +
    labs(
      title = "Volcano Plot",
      subtitle = paste0(contrast_name, " | red: ", group1_name, " | blue: ", group2_name, " | cut_off_pvalue = ", p_cut, " | cut_off_logFC = ", fc_cut),
      x = "log2(Fold Change)",
      y = "-log10(P-value)",
      colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

make_top_deg_barplot <- function(df, contrast_name, group1_name, group2_name) {
  plot_df <- df %>%
    filter(significant %in% c("Up", "Down")) %>%
    mutate(score = abs(log2FoldChange) * neg_log10_p) %>%
    arrange(desc(score)) %>%
    slice_head(n = 20) %>%
    arrange(log2FoldChange) %>%
    mutate(
      gene = factor(gene, levels = gene),
      display_group = ifelse(significant == "Up", group1_name, group2_name)
    )
  if (nrow(plot_df) == 0) {
    return(blank_plot("No significant genes", contrast_name))
  }
  ggplot(plot_df, aes(x = gene, y = log2FoldChange, fill = display_group)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = setNames(c("#ff4757", "#546de5"), c(group1_name, group2_name))) +
    labs(
      title = paste0("Top Differential Genes: ", contrast_name),
      x = NULL,
      y = "log2(Fold Change)",
      fill = "Higher in"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

make_fc_density_plot <- function(df, contrast_name, group1_name, group2_name) {
  plot_df <- add_group_display_columns(df, group1_name, group2_name)
  ggplot(plot_df, aes(x = log2FoldChange, fill = display_group, colour = display_group)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    scale_fill_manual(values = setNames(c("#9CC6DC", "#D9DEE3", "#F29AA3"), c(group2_name, "None", group1_name))) +
    scale_colour_manual(values = setNames(c("#546de5", "#AAB2BB", "#ff4757"), c(group2_name, "None", group1_name))) +
    labs(
      title = paste0("log2FC Distribution: ", contrast_name),
      x = "log2(Fold Change)",
      y = "Density",
      fill = "Higher in",
      colour = "Higher in"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

make_ma_plot <- function(df, contrast_name, group1_name, group2_name) {
  abundance_col <- get_abundance_col(df)
  if (is.null(abundance_col)) {
    return(blank_plot("Abundance information unavailable", contrast_name))
  }
  plot_df <- add_group_display_columns(df, group1_name, group2_name)
  if (abundance_col == "baseMean") {
    plot_df$x_value <- log10(plot_df[[abundance_col]] + 1)
    x_label <- "log10(baseMean + 1)"
  } else {
    plot_df$x_value <- plot_df[[abundance_col]]
    x_label <- "logCPM"
  }
  ggplot(plot_df, aes(x = x_value, y = log2FoldChange, colour = display_group)) +
    geom_point(alpha = 0.65, size = 2) +
    geom_hline(yintercept = c(-opt$cut_off_logFC, opt$cut_off_logFC), linetype = "dashed", linewidth = 0.45, colour = "#6B7280") +
    scale_colour_manual(values = setNames(c("#546de5", "#C7CDD4", "#ff4757"), c(group2_name, "None", group1_name))) +
    labs(
      title = paste0("MA-style Plot: ", contrast_name),
      x = x_label,
      y = "log2(Fold Change)",
      colour = "Higher in"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

convert_gene_ids <- function(genes, species_cfg) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  cache_key <- paste(class(species_cfg$orgdb), collapse = "_")
  if (!exists(cache_key, envir = annotation_key_cache, inherits = FALSE)) {
    valid_symbols <- keys(species_cfg$orgdb, keytype = "SYMBOL")
    assign(cache_key, unique(valid_symbols[!is.na(valid_symbols)]), envir = annotation_key_cache)
  }
  valid_symbols <- get(cache_key, envir = annotation_key_cache, inherits = FALSE)
  genes <- intersect(genes, valid_symbols)
  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  gene_map <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = species_cfg$orgdb
  )
  gene_map <- gene_map %>%
    filter(!is.na(SYMBOL), !is.na(ENTREZID)) %>%
    distinct(SYMBOL, ENTREZID)
  gene_map
}

plot_go_enrichment <- function(go_df, file_path, title_text) {
  if (nrow(go_df) == 0) {
    save_png_plot(blank_plot("No significant GO terms", title_text), file_path, width = 10, height = 7)
    return(invisible(NULL))
  }
  plot_df <- go_df %>%
    filter(!is.na(ONTOLOGY), !is.na(Description)) %>%
    arrange(p.adjust, pvalue) %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(
      Description = wrap_text(Description, 36),
      Description = factor(Description, levels = rev(unique(Description)))
    )
  p <- ggplot(plot_df, aes(x = Count, y = Description, size = Count, colour = p.adjust)) +
    geom_point(alpha = 0.9) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#D1495B", high = "#F6D6D9", trans = "reverse") +
    labs(
      title = title_text,
      x = "Gene Count",
      y = NULL,
      colour = "Adjusted P",
      size = "Gene Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text.y = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  save_png_plot(p, file_path, width = 10.5, height = 8)
}

plot_kegg_enrichment <- function(kegg_df, file_path, title_text) {
  if (nrow(kegg_df) == 0) {
    save_png_plot(blank_plot("No significant KEGG pathways", title_text), file_path, width = 10, height = 7)
    return(invisible(NULL))
  }
  plot_df <- kegg_df %>%
    filter(!is.na(Description)) %>%
    arrange(p.adjust, pvalue) %>%
    slice_head(n = 15) %>%
    mutate(
      Description = wrap_text(Description, 40),
      Description = factor(Description, levels = rev(unique(Description)))
    )
  p <- ggplot(plot_df, aes(x = Count, y = Description, size = Count, colour = p.adjust)) +
    geom_point(alpha = 0.92) +
    scale_colour_gradient(low = "#D1495B", high = "#F6D6D9", trans = "reverse") +
    labs(
      title = title_text,
      x = "Gene Count",
      y = NULL,
      colour = "Adjusted P",
      size = "Gene Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  save_png_plot(p, file_path, width = 10.5, height = 7.5)
}

run_ora <- function(gene_df, species_cfg, save_dir) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  gene_path <- file.path(save_dir, "diff_genes.csv")
  go_path <- file.path(save_dir, "GO_data.csv")
  kegg_path <- file.path(save_dir, "kegg_data.csv")
  write_csv_file(gene_df, gene_path)
  if (nrow(gene_df) == 0) {
    write_csv_file(data.frame(), go_path)
    write_csv_file(data.frame(), kegg_path)
    save_png_plot(blank_plot("No differential genes", basename(save_dir)), file.path(save_dir, "GO_enrich.png"), width = 10, height = 7)
    save_png_plot(blank_plot("No differential genes", basename(save_dir)), file.path(save_dir, "KEGG_enrich.png"), width = 10, height = 7)
    return(invisible(NULL))
  }
  gene_map <- convert_gene_ids(gene_df$gene, species_cfg)
  if (nrow(gene_map) == 0) {
    write_csv_file(data.frame(), go_path)
    write_csv_file(data.frame(), kegg_path)
    save_png_plot(blank_plot("Gene symbols not mapped", paste0("Check --spacies: ", opt$spacies)), file.path(save_dir, "GO_enrich.png"), width = 10, height = 7)
    save_png_plot(blank_plot("Gene symbols not mapped", paste0("Check --spacies: ", opt$spacies)), file.path(save_dir, "KEGG_enrich.png"), width = 10, height = 7)
    return(invisible(NULL))
  }
  go_result <- enrichGO(
    gene = unique(gene_map$ENTREZID),
    OrgDb = species_cfg$orgdb,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  kegg_result <- enrichKEGG(
    gene = unique(gene_map$ENTREZID),
    organism = species_cfg$kegg,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  go_df <- as.data.frame(go_result@result)
  kegg_df <- as.data.frame(kegg_result@result)
  write_csv_file(go_df, go_path)
  write_csv_file(kegg_df, kegg_path)
  plot_go_enrichment(go_df, file.path(save_dir, "GO_enrich.png"), paste0(basename(save_dir), " GO Enrichment"))
  plot_kegg_enrichment(kegg_df, file.path(save_dir, "KEGG_enrich.png"), paste0(basename(save_dir), " KEGG Enrichment"))
}

build_ranked_gene_list <- function(df, species_cfg) {
  ranked_df <- df %>%
    select(gene, log2FoldChange, pval_use) %>%
    distinct(gene, .keep_all = TRUE) %>%
    filter(!is.na(gene), gene != "", !is.na(log2FoldChange))
  gene_map <- convert_gene_ids(ranked_df$gene, species_cfg)
  if (nrow(gene_map) == 0) {
    return(numeric(0))
  }
  ranked_map <- ranked_df %>%
    inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
    mutate(rank_score = sign(log2FoldChange) * (-log10(ifelse(is.na(pval_use) | pval_use <= 0, 1e-300, pval_use)))) %>%
    group_by(ENTREZID) %>%
    summarise(rank_score = mean(rank_score), .groups = "drop") %>%
    filter(!is.na(rank_score), is.finite(rank_score)) %>%
    arrange(desc(rank_score))
  gene_list <- ranked_map$rank_score
  names(gene_list) <- ranked_map$ENTREZID
  sort(gene_list, decreasing = TRUE)
}

plot_gsea_summary <- function(gsea_df, file_path, title_text) {
  if (nrow(gsea_df) == 0) {
    save_png_plot(blank_plot("No enriched gene sets", title_text), file_path, width = 10, height = 7)
    return(invisible(NULL))
  }
  plot_df <- gsea_df %>%
    filter(!is.na(Description), !is.na(NES), !is.na(p.adjust)) %>%
    arrange(p.adjust, desc(abs(NES))) %>%
    slice_head(n = 15) %>%
    mutate(
      Description = wrap_text(Description, 38),
      Description = factor(Description, levels = rev(unique(Description))),
      direction = ifelse(NES >= 0, "Activated", "Suppressed")
    )
  p <- ggplot(plot_df, aes(x = NES, y = Description, size = setSize, colour = p.adjust, shape = direction)) +
    geom_point(alpha = 0.92) +
    scale_colour_gradient(low = "#D1495B", high = "#F6D6D9", trans = "reverse") +
    scale_shape_manual(values = c("Activated" = 16, "Suppressed" = 17)) +
    labs(
      title = title_text,
      x = "Normalized Enrichment Score",
      y = NULL,
      colour = "Adjusted P",
      size = "Set Size",
      shape = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  save_png_plot(p, file_path, width = 10.5, height = 7.5)
}

run_gsea <- function(df, species_cfg, save_dir) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  go_path <- file.path(save_dir, "GSEA_GO_data.csv")
  kegg_path <- file.path(save_dir, "GSEA_KEGG_data.csv")
  gene_list <- build_ranked_gene_list(df, species_cfg)
  if (length(gene_list) < 10 || length(unique(gene_list)) < 2) {
    write_csv_file(data.frame(), go_path)
    write_csv_file(data.frame(), kegg_path)
    save_png_plot(blank_plot("Insufficient ranked genes", "GSEA GO"), file.path(save_dir, "GSEA_GO_plot.png"), width = 10, height = 7)
    save_png_plot(blank_plot("Insufficient ranked genes", "GSEA KEGG"), file.path(save_dir, "GSEA_KEGG_plot.png"), width = 10, height = 7)
    return(invisible(NULL))
  }
  gsea_go <- gseGO(
    geneList = gene_list,
    OrgDb = species_cfg$orgdb,
    keyType = "ENTREZID",
    ont = "ALL",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  gsea_kegg <- gseKEGG(
    geneList = gene_list,
    organism = species_cfg$kegg,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  go_df <- as.data.frame(gsea_go@result)
  kegg_df <- as.data.frame(gsea_kegg@result)
  write_csv_file(go_df, go_path)
  write_csv_file(kegg_df, kegg_path)
  plot_gsea_summary(go_df, file.path(save_dir, "GSEA_GO_plot.png"), "GSEA GO")
  plot_gsea_summary(kegg_df, file.path(save_dir, "GSEA_KEGG_plot.png"), "GSEA KEGG")
}

annotate_gene_ownership <- function(df, group1_name, group2_name) {
  df$higher_in_group <- ifelse(
    df$log2FoldChange > 0,
    group1_name,
    ifelse(df$log2FoldChange < 0, group2_name, "None")
  )
  df$lower_in_group <- ifelse(
    df$log2FoldChange > 0,
    group2_name,
    ifelse(df$log2FoldChange < 0, group1_name, "None")
  )
  df
}

write_group_specific_outputs <- function(df, target_dir, species_cfg) {
  write_csv_file(df, file.path(target_dir, "diff_strict.csv"))
  write_csv_file(df, file.path(target_dir, "diff_loose.csv"))
  run_ora(df, species_cfg, target_dir)
}

save_contrast_heatmap <- function(df, file_path) {
  heat_df <- df %>%
    select(gene, contrast, log2FoldChange) %>%
    distinct(gene, contrast, .keep_all = TRUE)
  heat_wide <- heat_df %>%
    pivot_wider(names_from = contrast, values_from = log2FoldChange)
  if (nrow(heat_wide) < 2 || ncol(heat_wide) < 3) {
    save_png_plot(blank_plot("Heatmap unavailable", "Need at least two contrasts"), file_path, width = 10, height = 7)
    return(invisible(NULL))
  }
  heat_mat <- as.matrix(heat_wide[, -1, drop = FALSE])
  rownames(heat_mat) <- heat_wide$gene
  heat_mat[is.na(heat_mat)] <- 0
  max_abs <- apply(abs(heat_mat), 1, max)
  top_n <- min(60, length(max_abs))
  top_genes <- names(sort(max_abs, decreasing = TRUE))[seq_len(top_n)]
  heat_mat <- heat_mat[top_genes, , drop = FALSE]
  if (nrow(heat_mat) < 2 || ncol(heat_mat) < 2 || length(unique(as.vector(heat_mat))) < 2) {
    save_png_plot(blank_plot("Heatmap unavailable", "Insufficient variation"), file_path, width = 10, height = 7)
    return(invisible(NULL))
  }
  color_pal <- colorRampPalette(c("#2E86AB", "#F7F7F7", "#D1495B"))(100)
  pheatmap(
    mat = heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    color = color_pal,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize = 9,
    filename = file_path,
    width = 10.5,
    height = 8
  )
}

species_cfg <- get_species_config(opt$spacies)
aggregated_diffexp <- load_diff_results(opt$input_dir, opt$algorithm)
aggregated_diffexp <- ensure_gene_column(aggregated_diffexp)

contrast_list <- if ("contrast" %in% colnames(aggregated_diffexp)) unique(aggregated_diffexp$contrast) else "comparison"
summary_rows <- list()
is_marker_input <- identical(basename(opt$input_dir), "marker_genes_pval.csv")

for (contrast_name in contrast_list) {
  if ("contrast" %in% colnames(aggregated_diffexp)) {
    df_ct_raw <- aggregated_diffexp %>% filter(contrast == contrast_name)
  } else {
    df_ct_raw <- aggregated_diffexp
  }
  if (nrow(df_ct_raw) == 0) {
    next
  }
  df_ct <- prepare_diff_table(df_ct_raw, opt$cut_off_pvalue, opt$cut_off_logFC)
  contrast_groups <- resolve_group_names(df_ct_raw, contrast_name, opt$sample_list)
  group1_name <- contrast_groups$group1
  group2_name <- contrast_groups$group2
  df_ct <- annotate_gene_ownership(df_ct, group1_name, group2_name)
  ct_dir <- if (identical(contrast_name, contrast_list[1])) parent_dir else file.path(parent_dir, safe_name(contrast_name))
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv_file(df_ct, file.path(ct_dir, "diff_all.csv"))

  if (is_marker_input) {
    volcano_plot <- make_volcano_plot(df_ct, contrast_name, opt$cut_off_pvalue, opt$cut_off_logFC, group1_name, group2_name)
    save_png_plot(volcano_plot, file.path(ct_dir, "volcano.png"), width = 10.5, height = 8)
  }

  save_png_plot(make_top_deg_barplot(df_ct, contrast_name, group1_name, group2_name), file.path(ct_dir, "top_deg_barplot.png"), width = 10, height = 7.5)
  save_png_plot(make_fc_density_plot(df_ct, contrast_name, group1_name, group2_name), file.path(ct_dir, "log2fc_density.png"), width = 10, height = 7)
  save_png_plot(make_ma_plot(df_ct, contrast_name, group1_name, group2_name), file.path(ct_dir, "ma_plot.png"), width = 10, height = 7)

  positive_df <- df_ct %>% filter(significant == "Up")
  negative_df <- df_ct %>% filter(significant == "Down")
  positive_dir <- file.path(ct_dir, "positive")
  negative_dir <- file.path(ct_dir, "negative")
  positive_group_dir <- file.path(ct_dir, paste0("higher_in_", safe_name(group1_name)))
  negative_group_dir <- file.path(ct_dir, paste0("higher_in_", safe_name(group2_name)))
  gsea_dir <- file.path(ct_dir, "gsea")

  write_group_specific_outputs(positive_df, positive_dir, species_cfg)
  write_group_specific_outputs(negative_df, negative_dir, species_cfg)
  write_group_specific_outputs(positive_df, positive_group_dir, species_cfg)
  write_group_specific_outputs(negative_df, negative_group_dir, species_cfg)
  run_gsea(df_ct, species_cfg, gsea_dir)

  summary_rows[[paste0(contrast_name, "_", group1_name)]] <- data.frame(
    contrast = contrast_name,
    group_name = group1_name,
    direction = "Higher",
    count = nrow(positive_df),
    output_dir = basename(positive_group_dir)
  )
  summary_rows[[paste0(contrast_name, "_", group2_name)]] <- data.frame(
    contrast = contrast_name,
    group_name = group2_name,
    direction = "Higher",
    count = nrow(negative_df),
    output_dir = basename(negative_group_dir)
  )
}

summary_df <- bind_rows(summary_rows)
write_csv_file(summary_df, file.path(parent_dir, "contrast_summary.csv"))

if (nrow(summary_df) > 0) {
  summary_plot <- ggplot(summary_df, aes(x = contrast, y = count, fill = group_name)) +
    geom_col(position = "dodge", width = 0.72) +
    labs(
      title = "Differential Expression Summary",
      x = "Contrast",
      y = "Gene Count",
      fill = "Higher expression in"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 28, hjust = 1),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  save_png_plot(summary_plot, file.path(parent_dir, "contrast_summary.png"), width = 11, height = 7)
}

if ("contrast" %in% colnames(aggregated_diffexp)) {
  prepared_all <- aggregated_diffexp %>%
    group_by(contrast) %>%
    group_modify(~prepare_diff_table(.x, opt$cut_off_pvalue, opt$cut_off_logFC)) %>%
    ungroup()
  save_contrast_heatmap(prepared_all, file.path(parent_dir, "contrast_log2fc_heatmap.png"))
}

if (!file.exists(opt$output_path)) {
  write_csv_file(data.frame(), opt$output_path)
}
