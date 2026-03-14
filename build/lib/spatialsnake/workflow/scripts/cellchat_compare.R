library(optparse)
library(CellChat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)

option_list <- list(
  make_option(c("--input_rds1"), type="character"),
  make_option(c("--input_rds2"), type="character", default=""),
  make_option(c("--sample_name1"), type="character", default="group1"),
  make_option(c("--sample_name2"), type="character", default="group2"),
  make_option(c("-o", "--output_dir"), type="character"),
  make_option(c("--pathways"), type="character", default=""),
  make_option(c("--source_cells"), type="character", default=""),
  make_option(c("--target_cells"), type="character", default=""),
  make_option(c("--receiver_cells"), type="character", default=""),
  make_option(c("--bubble_angle"), type="double", default=45),
  make_option(c("--bubble_remove_isolate"), type="logical", default=TRUE),
  make_option(c("--do_ranknet"), type="logical", default=TRUE),
  make_option(c("--do_role_heatmap"), type="logical", default=TRUE),
  make_option(c("--do_pathway_plots"), type="logical", default=TRUE),
  make_option(c("--do_compare_overview"), type="logical", default=TRUE),
  make_option(c("--do_compare_bubble"), type="logical", default=TRUE),
  make_option(c("--do_single_bubble"), type="logical", default=TRUE),
  make_option(c("--do_gene_expression"), type="logical", default=FALSE),
  make_option(c("--gene_colors"), type="character", default="white,#FEC44F,#D95F0E"),
  make_option(c("--gene_plot_type"), type="character", default="dot"),
  make_option(c("--save_merged"), type="logical", default=TRUE)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds1) || opt$input_rds1 == "") {
  stop("input_rds1 is required")
}
if (is.null(opt$output_dir) || opt$output_dir == "") {
  stop("output_dir is required")
}

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

trim_vec <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

split_vec <- function(x) {
  if (is.null(x) || x == "" || tolower(x) %in% c("null", "none")) {
    return(NULL)
  }
  trim_vec(unlist(strsplit(x, ",")))
}

parse_cells <- function(value, levels_vec) {
  vals <- split_vec(value)
  if (is.null(vals) || length(vals) == 0) {
    return(NULL)
  }
  if (all(grepl("^[0-9]+$", vals))) {
    return(as.integer(vals))
  }
  idx <- match(vals, levels_vec)
  if (any(is.na(idx))) {
    stop(paste0("Unknown cell types: ", paste(vals[is.na(idx)], collapse = ",")))
  }
  idx
}

resolve_name <- function(path, fallback) {
  if (!is.null(path) && path != "" && !is.null(fallback) && fallback != "") {
    return(fallback)
  }
  tools::file_path_sans_ext(basename(path))
}

format_vec <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("无")
  }
  paste(x, collapse = ", ")
}

extract_lr_pairs <- function(cellchat_obj) {
  net_lr <- subsetCommunication(cellchat_obj)
  if (is.null(net_lr) || nrow(net_lr) == 0) {
    return(NULL)
  }
  if (all(c("ligand", "receptor") %in% names(net_lr))) {
    return(unique(paste(net_lr$ligand, net_lr$receptor, sep = "_")))
  }
  if ("interaction_name" %in% names(net_lr)) {
    return(unique(net_lr$interaction_name))
  }
  if ("interaction_name_2" %in% names(net_lr)) {
    return(unique(net_lr$interaction_name_2))
  }
  NULL
}

print_visualization_options <- function(title, cell_levels, pathway_list, lr_pairs) {
  message("=== 可视化可选项: ", title, " ===")
  message("细胞类型(", length(cell_levels), "): ", format_vec(cell_levels))
  message("通路(", length(pathway_list), "): ", format_vec(pathway_list))
  message("配受体(", length(lr_pairs), "): ", format_vec(lr_pairs))
}

plot_gene_expression <- function(cellchat, pathways, file_prefix, color_vec, plot_type) {
  if (is.null(pathways) || length(pathways) == 0) {
    return()
  }
  pdf(file.path(opt$output_dir, paste0(file_prefix, "_gene_expression.pdf")), width = 10, height = 6)
  plotGeneExpression(cellchat, signaling = pathways, type = plot_type, col = color_vec)
  dev.off()
}

plot_pathway_outputs <- function(cellchat, pathways, receiver_cells, file_prefix) {
  if (is.null(pathways) || length(pathways) == 0) {
    return()
  }
  pdf(file.path(opt$output_dir, paste0(file_prefix, "_aggregate_hierarchy.pdf")), width = 10, height = 8)
  netVisual_aggregate(cellchat, signaling = pathways, vertex.receiver = receiver_cells, layout = "hierarchy")
  dev.off()
  pdf(file.path(opt$output_dir, paste0(file_prefix, "_aggregate_chord.pdf")), width = 10, height = 8)
  netVisual_aggregate(cellchat, signaling = pathways, vertex.receiver = receiver_cells, layout = "chord")
  dev.off()
  pdf(file.path(opt$output_dir, paste0(file_prefix, "_pathway_heatmap.pdf")), width = 10, height = 8)
  netVisual_heatmap(cellchat, signaling = pathways, color.heatmap = "Reds")
  dev.off()
  pdf(file.path(opt$output_dir, paste0(file_prefix, "_pathway_contribution.pdf")), width = 10, height = 8)
  for (pathway in pathways) {
    netAnalysis_contribution(cellchat, signaling = pathway)
  }
  dev.off()
}

plot_single_bubble <- function(cellchat, source_cells, target_cells, pathways, file_prefix) {
  if (is.null(source_cells)) {
    return()
  }
  pdf(file.path(opt$output_dir, paste0(file_prefix, "_bubble.pdf")), width = 12, height = 6)
  if (is.null(target_cells)) {
    netVisual_bubble(cellchat, sources.use = source_cells, signaling = pathways, remove.isolate = opt$bubble_remove_isolate)
  } else {
    netVisual_bubble(cellchat, sources.use = source_cells, targets.use = target_cells, signaling = pathways, remove.isolate = opt$bubble_remove_isolate)
  }
  dev.off()
}

plot_compare_bubble <- function(cellchat_merged, source_cells, target_cells) {
  if (is.null(source_cells) || is.null(target_cells)) {
    return()
  }
  p1 <- netVisual_bubble(cellchat_merged, sources.use = source_cells, targets.use = target_cells, comparison = c(1, 2), max.dataset = 2, angle.x = opt$bubble_angle, remove.isolate = opt$bubble_remove_isolate)
  p2 <- netVisual_bubble(cellchat_merged, sources.use = source_cells, targets.use = target_cells, comparison = c(1, 2), max.dataset = 1, angle.x = opt$bubble_angle, remove.isolate = opt$bubble_remove_isolate)
  pc <- p1 + p2
  ggsave(file.path(opt$output_dir, "compare_lr_regulated.pdf"), pc, width = 12, height = 5.5)
}

pathways <- split_vec(opt$pathways)
gene_colors <- split_vec(opt$gene_colors)

if (!is.null(opt$input_rds2) && opt$input_rds2 != "") {
  cellchat1 <- readRDS(opt$input_rds1)
  cellchat2 <- readRDS(opt$input_rds2)
  name1 <- resolve_name(opt$input_rds1, opt$sample_name1)
  name2 <- resolve_name(opt$input_rds2, opt$sample_name2)
  if (length(levels(cellchat1@idents))!=length(levels(cellchat2@idents))) {
    if (length(levels(cellchat1@idents)) < length(levels(cellchat2@idents))) {
      group.new = levels(cellchat2@idents)
      cellchat1 <- liftCellChat(cellchat1, group.new)
    } else {
      group.new = levels(cellchat1@idents)
      cellchat2 <- liftCellChat(cellchat2, group.new)
    }
  }
  
  cco.list <- list()
  cco.list[[name1]] <- cellchat1
  cco.list[[name2]] <- cellchat2
  cellchat_merged <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
  pathway_union <- union(cco.list[[name1]]@netP$pathways, cco.list[[name2]]@netP$pathways)
  lr_pairs_union <- unique(c(extract_lr_pairs(cco.list[[name1]]), extract_lr_pairs(cco.list[[name2]])))
  print_visualization_options(paste0(name1, " vs ", name2), levels(cellchat_merged@idents), pathway_union, lr_pairs_union)
  if (opt$save_merged) {
    saveRDS(cco.list, file.path(opt$output_dir, "cellchat_object.list.rds"))
    saveRDS(cellchat_merged, file.path(opt$output_dir, "cellchat_merged.rds"))
  }
  if (opt$do_compare_overview) {
    gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "count")
    gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
    p <- gg1 + gg2
    ggsave(file.path(opt$output_dir, "compare_overview_number_strength.pdf"), p, width = 6, height = 4)
    pdf(file.path(opt$output_dir, "compare_diff_number_strength_net.pdf"), width = 10, height = 5)
    par(mfrow = c(1, 2))
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
    dev.off()
    pdf(file.path(opt$output_dir, "compare_heatmap_count_weight.pdf"), width = 10, height = 6)
    h1 <- netVisual_heatmap(cellchat_merged)
    h2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
    print(h1 + h2)
    dev.off()
  }
  if (opt$do_ranknet) {
    gg1 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
    p <- gg1 + gg2
    ggsave(file.path(opt$output_dir, "compare_pathway_strength.pdf"), p, width = 10, height = 6)
  }
  if (opt$do_role_heatmap) {
    cco.list[[name1]] <- netAnalysis_computeCentrality(cco.list[[name1]], slot.name = "netP")
    cco.list[[name2]] <- netAnalysis_computeCentrality(cco.list[[name2]], slot.name = "netP")
    for (pattern in c("all", "outgoing", "incoming")) {
      ht1 <- netAnalysis_signalingRole_heatmap(cco.list[[name1]], pattern = pattern, signaling = pathway_union, title = name1, width = 8, height = 10)
      ht2 <- netAnalysis_signalingRole_heatmap(cco.list[[name2]], pattern = pattern, signaling = pathway_union, title = name2, width = 8, height = 10)
      pdf(file.path(opt$output_dir, paste0("compare_signaling_role_", pattern, ".pdf")), width = 16, height = 10)
      draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
      dev.off()
    }
  }
  if (opt$do_compare_bubble) {
    source_cells <- parse_cells(opt$source_cells, levels(cellchat_merged@idents))
    target_cells <- parse_cells(opt$target_cells, levels(cellchat_merged@idents))
    plot_compare_bubble(cellchat_merged, source_cells, target_cells)
  }
  if (opt$do_pathway_plots) {
    receiver_cells <- parse_cells(opt$receiver_cells, levels(cellchat_merged@idents))
    plot_pathway_outputs(cellchat_merged, pathways, receiver_cells, "compare_pathway")
  }
  if (opt$do_gene_expression) {
    plot_gene_expression(cellchat_merged, pathways, "compare_pathway", gene_colors, opt$gene_plot_type)
  }
} else {
  cellchat <- readRDS(opt$input_rds1)
  sample_name <- resolve_name(opt$input_rds1, opt$sample_name1)
  receiver_cells <- parse_cells(opt$receiver_cells, levels(cellchat@idents))
  source_cells <- parse_cells(opt$source_cells, levels(cellchat@idents))
  target_cells <- parse_cells(opt$target_cells, levels(cellchat@idents))
  lr_pairs_single <- extract_lr_pairs(cellchat)
  print_visualization_options(sample_name, levels(cellchat@idents), cellchat@netP$pathways, lr_pairs_single)
  if (opt$do_pathway_plots) {
    plot_pathway_outputs(cellchat, pathways, receiver_cells, sample_name)
  }
  if (opt$do_single_bubble) {
    plot_single_bubble(cellchat, source_cells, target_cells, pathways, sample_name)
  }
  if (opt$do_gene_expression) {
    plot_gene_expression(cellchat, pathways, sample_name, gene_colors, opt$gene_plot_type)
  }
}
