library(optparse)
library(Seurat)
library(spacexr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)

option_list <- list(
  make_option(c("-s", "--spatial_input"), type="character", help="Path to spatial data (Seurat .rds or .h5Seurat or .h5ad)"),
  make_option(c("-r", "--sc_input"), type="character", help="Path to single cell data (Seurat .rds or .h5Seurat or .h5ad)"),
  make_option(c("-o", "--output_dir"), type="character", help="Directory to save results"),
  make_option(c("--sample_id"), type="character", default="sample", help="Sample ID"),
  make_option(c("--mode"), type="character", default="doublet", help="RCTD mode: 'doublet' or 'full'"),
  make_option(c("--max_cores"), type="integer", default=64, help="Number of cores for RCTD"),
  make_option(c("--sc_cell_type_col"), type="character", default="celltype", help="Column name for cell types in SC data"),
  make_option(c("--spatial_cell_type_col"), type="character", default="celltype", help="Column name for cell types or clusters in spatial metadata"),
  make_option(c("--group_by"), type="character", default="sample", help="Column in Spatial Seurat object to group by (X-axis) for heatmap")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt$max_cores)

if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}


h5ad_loading <- function(input_h5ad, is_spatial = TRUE) {
  obj <- tryCatch({
    message("Attempting standard spatial conversion...")
    schard::h5ad2seurat_spatial(input_h5ad, use.raw = TRUE)
  }, error = function(e) {
    message("Standard spatial conversion failed or encountered error.")
    message(paste("Error message:", e$message))
    message("Attempting multi-sample conversion (simplify = FALSE) with merge loop...")
    
    visl <- schard::h5ad2seurat_spatial(input_h5ad, use.raw = TRUE, simplify = FALSE)
    
    if (length(visl) == 0) {
      stop("Multi-sample conversion returned no objects.")
    }
    
    left0 <- NULL
    for (i in visl) {
      if (is.null(left0)) {
         left0 <- merge(i, left0) 
      } else {
         left0 <- merge(i, left0)
      }
    }
    return(left0)
  })
  return(obj)
}






load_seurat <- function(file_path, is_spatial = FALSE) {
  if (is_spatial) {
    if (grepl("\\.h5ad$", file_path, ignore.case = TRUE)) {
      return(h5ad_loading(file_path, is_spatial = TRUE))
    } else {
      return(readRDS(file_path))
    }
  } else {
    if (grepl("\\.rds$", file_path, ignore.case = TRUE)) {
      return(readRDS(file_path))
    } else if (grepl("\\.h5ad$", file_path, ignore.case = TRUE)) {
      message("Detected h5ad format for Single Cell Data. Converting to Seurat using schard...")
      if (!requireNamespace("schard", quietly = TRUE)) stop("schard package required for reading h5ad")
      return(schard::h5ad2seurat(file_path, use.raw = TRUE))
    } else {
      stop("Unsupported file format for Single Cell Data: ", file_path)
    }
  }
}

plot_rctd_dotplot <- function(weights, spatial_meta, unsupervised_cluster_col, save_path, sample_id) {
  if (!unsupervised_cluster_col %in% colnames(spatial_meta)) {
    stop(sprintf("Column '%s' not found in spatial metadata", unsupervised_cluster_col))
  }

  common_spots <- intersect(rownames(spatial_meta), rownames(weights))
  if (length(common_spots) < 3) {
    stop("Too few overlapping spots between spatial metadata and RCTD full weights.")
  }

  spatial_meta <- spatial_meta[common_spots, , drop = FALSE]
  weights <- weights[common_spots, , drop = FALSE]

  cluster_vec <- spatial_meta[[unsupervised_cluster_col]]
  valid_spots <- !is.na(cluster_vec) & as.character(cluster_vec) != ""
  cluster_vec <- as.character(cluster_vec[valid_spots])
  weights <- weights[valid_spots, , drop = FALSE]

  if (nrow(weights) < 3) {
    stop("Too few valid spots remain after filtering NA cluster labels.")
  }

  cluster_factor <- factor(cluster_vec)
  if (length(levels(cluster_factor)) < 2) {
    stop("At least two cluster groups are required to compute the full-mode correlation dotplot.")
  }

  cluster_design <- stats::model.matrix(~ 0 + cluster_factor)
  colnames(cluster_design) <- sub("^cluster_factor", "", colnames(cluster_design))

  weights_dense <- as.matrix(weights)
  mode(weights_dense) <- "numeric"

  keep_celltypes <- apply(weights_dense, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  weights_dense <- weights_dense[, keep_celltypes, drop = FALSE]

  if (ncol(weights_dense) == 0) {
    stop("All full-mode RCTD weight columns have zero variance.")
  }

  cor_mat <- stats::cor(
    weights_dense,
    cluster_design,
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("cell_type", "cluster", "correlation")
  cor_df$abs_correlation <- abs(cor_df$correlation)

  celltype_order <- names(sort(tapply(cor_df$abs_correlation, cor_df$cell_type, max), decreasing = TRUE))
  cluster_order <- names(sort(tapply(cor_df$abs_correlation, cor_df$cluster, max), decreasing = TRUE))

  cor_df$cell_type <- factor(cor_df$cell_type, levels = rev(celltype_order))
  cor_df$cluster <- factor(cor_df$cluster, levels = cluster_order)
  cor_df$label <- ifelse(cor_df$abs_correlation >= 0.30, sprintf("%.2f", cor_df$correlation), "")

  fig_width <- max(8, length(cluster_order) * 0.7 + 3)
  fig_height <- max(8, length(celltype_order) * 0.28 + 3)

  p <- ggplot(cor_df, aes(x = .data$cluster, y = .data$cell_type)) +
    geom_point(
      aes(size = .data$abs_correlation, fill = .data$correlation),
      shape = 21,
      colour = "grey20",
      stroke = 0.2,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = .data$label),
      size = 2.6,
      fontface = "bold",
      colour = "black"
    ) +
    scale_size_continuous(
      range = c(1.8, 9.5),
      limits = c(0, 1),
      breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
      name = "|Pearson r|"
    ) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "#F7F7F7",
      high = "#B40426",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.5),
      name = "Pearson r"
    ) +
    labs(
      title = "Correlation between unsupervised clusters and RCTD full deconvolution",
      subtitle = sprintf(
        "%s | RCTD full mode | matrix = normalized full cell proportion weights",
        sample_id
      ),
      x = unsupervised_cluster_col,
      y = "RCTD cell type",
      caption = "Bubble size indicates absolute correlation; color indicates signed correlation."
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      plot.title = element_text(face = "bold", size = 16, colour = "black"),
      plot.subtitle = element_text(size = 11, colour = "#4D4D4D"),
      plot.caption = element_text(size = 10, colour = "#4D4D4D"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      legend.position = "right",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "#FBFBFB", colour = NA)
    )

  ggsave(
    filename = save_path,
    plot = p,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 400,
    bg = "white"
  )

  invisible(cor_mat)
}

message("Loading Single Cell Data from: ", opt$sc_input)
merge_seurat_object <- load_seurat(opt$sc_input, is_spatial = FALSE)

message("Loading Spatial Data from: ", opt$spatial_input)
left0 <- load_seurat(opt$spatial_input, is_spatial = TRUE)

message("Preparing Reference...")
sc_cell_type_col <- opt$sc_cell_type_col
spatial_cell_type_col <- opt$spatial_cell_type_col
if (!sc_cell_type_col %in% colnames(merge_seurat_object@meta.data)) {
  stop(paste("Single-cell cell type column", sc_cell_type_col, "not found in SC metadata"))
}

celltype_all <- as.factor(merge_seurat_object@meta.data[[sc_cell_type_col]])
names(celltype_all) <- colnames(merge_seurat_object)
valid_cells <- !is.na(celltype_all) & as.character(celltype_all) != ""
if (!all(valid_cells)) {
  merge_seurat_object <- subset(merge_seurat_object, cells = names(celltype_all)[valid_cells])
  celltype_all <- as.factor(merge_seurat_object@meta.data[[sc_cell_type_col]])
  names(celltype_all) <- colnames(merge_seurat_object)
}
celltype_count <- sort(table(celltype_all), decreasing = TRUE)
keep_celltypes <- names(celltype_count[celltype_count >= 25])
if (length(keep_celltypes) == 0) {
  stop("No cell types with at least 25 cells found in the single-cell reference.")
}
merge_seurat_object <- subset(
  merge_seurat_object,
  cells = colnames(merge_seurat_object)[merge_seurat_object@meta.data[[sc_cell_type_col]] %in% keep_celltypes]
)
message("Retained cell types (>=25 cells): ", paste(keep_celltypes, collapse = ", "))
message("Retained reference cells: ", ncol(merge_seurat_object))

if (spatial_cell_type_col %in% colnames(left0@meta.data)) {
  Idents(left0) <- left0@meta.data[[spatial_cell_type_col]]
}

celltype <- as.factor(merge_seurat_object@meta.data[[sc_cell_type_col]])
names(celltype) <- colnames(merge_seurat_object)

if ("RNA" %in% names(merge_seurat_object@assays)) {
  assay_obj <- merge_seurat_object[["RNA"]]
  if (inherits(assay_obj, "Assay5")) {
    sc_counts <- assay_obj@layers$counts
    rownames(sc_counts) <- rownames(merge_seurat_object)
    colnames(sc_counts) <- colnames(merge_seurat_object)
  } else {
    sc_counts <- GetAssayData(merge_seurat_object, assay = "RNA", slot = "counts")
  }
} else {
  # 尝试使用默认 assay
  sc_counts <- GetAssayData(merge_seurat_object, slot = "counts")
}

celltype <- celltype[colnames(sc_counts)]

reference <- spacexr::Reference(sc_counts, celltype)

# 4. 准备 Spatial Data (空间数据)
message("Preparing Spatial Data...")
# 确保使用 Spatial assay (如果存在)
if ("Spatial" %in% names(left0@assays)) {
  DefaultAssay(left0) <- "Spatial"
}

# 获取 Spatial counts
if (inherits(left0[[DefaultAssay(left0)]], "Assay5")) {
  # Seurat v5: 尝试 JoinLayers 或者直接获取 counts
  st_count <- tryCatch({
    Seurat::JoinLayers(left0)
    Seurat::LayerData(left0, layer = "counts")
  }, error = function(e) {
    # 如果失败，尝试直接获取
    left0[[DefaultAssay(left0)]]$counts
  })
} else {
  st_count <- GetAssayData(left0, slot = "counts")
}


if (is.null(rownames(st_count))) rownames(st_count) <- rownames(left0)
if (is.null(colnames(st_count))) colnames(st_count) <- colnames(left0)


img_names <- Images(left0)
img_name <- if (length(img_names) > 0) img_names[1] else NULL
if (!is.null(img_name)) {
  message("Automatically detected image: ", img_name)
}

if ("boundaries" %in% slotNames(left0) && !is.null(tryCatch(left0@boundaries$centroids, error=function(e) NULL))) {
    st_loc <- left0@boundaries$centroids@coords
  } else {
    st_loc <- Seurat::GetTissueCoordinates(left0, image = img_name)
  }

if (is.null(st_loc)) {
  if (all(c("x", "y") %in% colnames(left0@meta.data))) {
    st_loc <- left0@meta.data[, c("x", "y")]
  } else if ("spatial" %in% names(left0@reductions)) {
    st_loc <- Embeddings(left0, "spatial")
  }
}

if (is.null(st_loc)) {
  stop("Could not find spatial coordinates in Seurat object.")
}

st_loc <- as.data.frame(st_loc)
colnames(st_loc) <- c("x", "y") # spacexr 要求列名为 x, y
rownames(st_loc) <- colnames(left0) # 确保行名是 spot/cell ID

# 创建 SpatialRNA 对象
query <- spacexr::SpatialRNA(st_loc, st_count)

# 5. 运行 RCTD
message("Running RCTD (mode: ", opt$mode, ")...")
RCTD <- create.RCTD(query, reference, max_cores = opt$max_cores)
RCTD <- run.RCTD(RCTD, doublet_mode = opt$mode)

# 6. 保存结果
message("Saving results to ", opt$output_dir)

# 保存完整结果对象 (可选，如果文件太大可以注释掉)
# saveRDS(RCTD, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_object.rds")))

# 提取结果 DataFrame
results_df <- RCTD@results$results_df
# 归一化权重 (Proportions)
norm_weights <- normalize_weights(RCTD@results$weights)

# 写入 CSV (包含 spot 列名)
write.csv(results_df, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_results.csv")), row.names = TRUE)
write.csv(norm_weights, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_weights.csv")), row.names = TRUE)


message("Generating plots...")
if (!is.null(results_df)) {
  left0 <- AddMetaData(left0, metadata = results_df)
}

common_weight_spots <- intersect(colnames(left0), rownames(norm_weights))
left0$RCTD_max_celltype <- NA_character_
left0$RCTD_max_weight <- NA_real_
if (length(common_weight_spots) > 0) {
  aligned_weights <- as.matrix(norm_weights[common_weight_spots, , drop = FALSE])
  max_idx <- max.col(aligned_weights, ties.method = "first")
  left0$RCTD_max_celltype[common_weight_spots] <- colnames(aligned_weights)[max_idx]
  left0$RCTD_max_weight[common_weight_spots] <- apply(aligned_weights, 1, max)
}

if ("first_type" %in% colnames(left0@meta.data)) {
  left0$first_type <- as.character(left0$first_type)
  left0$first_type[is.na(left0$first_type)] <- "Unknown"
}

saveRDS(left0, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_seurat.rds")))
saveRDS(RCTD, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD.rds")))
mode_label <- tolower(opt$mode)
min_count_to_show_label <- 20

if (identical(mode_label, "full")) {
  dotplot_col <- spatial_cell_type_col
  if (!dotplot_col %in% colnames(left0@meta.data) && opt$group_by %in% colnames(left0@meta.data)) {
    warning(
      paste0(
        "Column '", dotplot_col, "' not found in spatial metadata. Falling back to '",
        opt$group_by, "' for full-mode dotplot."
      )
    )
    dotplot_col <- opt$group_by
  }

  if (dotplot_col %in% colnames(left0@meta.data)) {
    dotplot_file <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_full_dotplot.png"))
    plot_rctd_dotplot(
        weights = norm_weights,
        spatial_meta = left0@meta.data,
        unsupervised_cluster_col = dotplot_col,
        save_path = dotplot_file,
        sample_id = opt$sample_id
      )
  } else {
    warning(
      paste0(
        "Neither '", spatial_cell_type_col, "' nor fallback column '", opt$group_by,
        "' was found in spatial metadata. Skipping full-mode dotplot."
      )
    )
    message("Available columns: ", paste(colnames(left0@meta.data), collapse = ", "))
  }
} else {
  if ("first_type" %in% colnames(left0@meta.data) && opt$group_by %in% colnames(left0@meta.data)) {
    message("Generating per-sample distribution plots...")

    sample_cell_counts <- left0@meta.data %>%
      dplyr::group_by(.data[[opt$group_by]], first_type) %>%
      dplyr::summarise(n = n(), .groups = "drop")

    p_sample <- ggplot(sample_cell_counts, aes(x = .data[[opt$group_by]], y = n, fill = first_type)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(
        aes(label = ifelse(n >= min_count_to_show_label, as.character(first_type), "")),
        position = position_stack(vjust = 0.5),
        size = 2
      ) +
      xlab(opt$group_by) +
      ylab("# of Spots") +
      ggtitle(paste0("Distribution of Cell Types across ", opt$group_by)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(
      filename = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_sample_dist_plot.png")),
      plot = p_sample,
      width = 12,
      height = 8
    )
  }

  group_col <- opt$group_by
  plot_cell_types <- function(data, label) {
    ggplot(data, aes(x = .data[[label]], y = n, fill = first_type)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(
        aes(label = ifelse(n >= min_count_to_show_label, as.character(first_type), "")),
        position = position_stack(vjust = 0.5),
        size = 2
      ) +
      xlab(label) +
      ylab("# of Spots") +
      ggtitle(paste0("Distribution of Cell Types across ", label)) +
      theme_minimal()
  }

  if ("first_type" %in% colnames(left0@meta.data) && group_col %in% colnames(left0@meta.data)) {
    message(paste0("Generating heatmap grouped by ", group_col, "..."))
    cell_type_banksy_counts <- left0[[]] %>%
      dplyr::count(first_type, .data[[group_col]], name = "n")
    colnames(cell_type_banksy_counts)[2] <- group_col
    p <- plot_cell_types(cell_type_banksy_counts, group_col)
    ggsave(
      filename = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_cluster_plot.png")),
      plot = p,
      width = 10,
      height = 7,
      units = "in"
    )

    counts <- table(left0$first_type, left0@meta.data[[group_col]])
    props <- prop.table(counts, margin = 2) * 100

    heatmap_file_png <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_heatmap.png"))
    heatmap_file_pdf <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_heatmap.pdf"))

    
    
    annotation_col <- data.frame(Group = colnames(props))
    rownames(annotation_col) <- colnames(props)
    names(annotation_col) <- group_col

    annotation_row <- data.frame(RCTD_Type = rownames(props))
    rownames(annotation_row) <- rownames(props)

    get_colors <- function(n) {
      if (n <= 9) return(RColorBrewer::brewer.pal(max(3, n), "Set1")[1:n])
      if (n <= 12) return(RColorBrewer::brewer.pal(n, "Set3"))
      rainbow(n)
    }

    anno_colors <- list()
    anno_colors[[group_col]] <- setNames(get_colors(ncol(props)), colnames(props))
    anno_colors[["RCTD_Type"]] <- setNames(get_colors(nrow(props)), rownames(props))

    tryCatch({
      png(filename = heatmap_file_png, width = 3200, height = 2400, res = 300)
      pheatmap(
        props,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        display_numbers = FALSE,
        color = reds_palette,
        border_color = "white",
        angle_col = 45,
        fontsize = 10,
        fontsize_row = 10,
        fontsize_col = 10,
        treeheight_row = 30,
        treeheight_col = 30,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = anno_colors,
        name = "Proportion (%)",
        main = paste0("RCTD first_type proportion by ", group_col)
      )
      dev.off()

      pdf(file = heatmap_file_pdf, width = 10, height = 8)
      pheatmap(
        props,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        display_numbers = FALSE,
        color = reds_palette,
        border_color = "white",
        angle_col = 45,
        fontsize = 10,
        fontsize_row = 10,
        fontsize_col = 10,
        treeheight_row = 30,
        treeheight_col = 30,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = anno_colors,
        name = "Proportion (%)",
        main = paste0("RCTD first_type proportion by ", group_col)
      )
      dev.off()

      message(paste0("Heatmap saved to ", heatmap_file_png))
    }, error = function(e) {
      message("Error generating heatmap: ", e$message)
      try({dev.off()}, silent = TRUE)
    })
  } else {
    warning(paste0("Column '", group_col, "' or 'first_type' not found in spatial metadata. Skipping heatmap."))
    message("Available columns: ", paste(colnames(left0@meta.data), collapse = ", "))
  }

  if (all(c("spot_class", "first_type") %in% colnames(left0@meta.data))) {
    spot_class_counts <- left0@meta.data %>%
      dplyr::select(first_type, spot_class) %>%
      dplyr::group_by(first_type, spot_class) %>%
      dplyr::summarise(count = n(), .groups = "drop") %>%
      dplyr::group_by(first_type) %>%
      dplyr::mutate(proportion = count / sum(count) * 100)

    spot_class_colors <- c(
      "singlet" = "#377EB8",
      "doublet_certain" = "#FF7F00",
      "doublet_uncertain" = "#984EA3",
      "reject" = "#E41A1C"
    )

    p_bar <- ggplot(spot_class_counts, aes(x = first_type, y = proportion, fill = spot_class)) +
      geom_bar(stat = "identity", position = "stack", width = 0.8) +
      scale_fill_manual(values = spot_class_colors) +
      labs(
        title = "Spot Class Distribution by Predicted Cell Type",
        x = "Predicted Cell Type (first_type)",
        y = "Proportion (%)",
        fill = "Spot Class"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )

    bar_plot_file <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_spot_class_bar.png"))
    ggsave(filename = bar_plot_file, plot = p_bar, width = 12, height = 6, dpi = 300)
    message(paste0("Spot class bar plot saved to ", bar_plot_file))
  }
}
message("RCTD Analysis Completed.")
