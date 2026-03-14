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
  make_option(c("--cell_type_col"), type="character", default="celltype", help="Column name for cell types in SC data"),
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

message("Loading Single Cell Data from: ", opt$sc_input)
merge_seurat_object <- load_seurat(opt$sc_input, is_spatial = FALSE)

message("Loading Spatial Data from: ", opt$spatial_input)
left0 <- load_seurat(opt$spatial_input, is_spatial = TRUE)

message("Preparing Reference...")
cell_type_col <- opt$cell_type_col
if (!cell_type_col %in% colnames(merge_seurat_object@meta.data)) {
  stop(paste("Cell type column", cell_type_col, "not found in SC metadata"))
}

Idents(left0) <- left0@meta.data[["celltype"]]
left0 <- subset(left0, downsample = 40000)

celltype <- as.factor(merge_seurat_object@meta.data[[cell_type_col]])
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
left0 <- AddMetaData(left0, metadata = results_df)
left0$first_type <- as.character(left0$first_type)
left0$first_type[is.na(left0$first_type)] <- "Unknown"
saveRDS(left0, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_seurat.rds")))


min_count_to_show_label <- 20

# New: Plot per-sample distribution if 'sample' column exists
if (opt$group_by %in% colnames(left0@meta.data)) {
  message("Generating per-sample distribution plots...")

  sample_cell_counts <- left0@meta.data %>%
    dplyr::group_by(sample, first_type) %>%
    dplyr::summarise(n = n(), .groups = 'drop')

  p_sample <- ggplot(sample_cell_counts, aes(x = sample, y = n, fill = first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, as.character(first_type), "")), 
              position = position_stack(vjust = 0.5), size = 2) +
    xlab("Sample") +
    ylab("# of Spots") +
    ggtitle("Distribution of Cell Types across Samples") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  ggsave(
    filename = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_sample_dist_plot.png")),
    plot = p_sample,
    width = 12,   
    height = 8
  )
}

group_col <- opt$cell_type_col
plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, as.character(first_type), "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}
if (group_col %in% colnames(left0@meta.data)) {
  message(paste0("Generating heatmap grouped by ", group_col, "..."))
  cell_type_banksy_counts <- left0[[]] %>%
  dplyr::count(first_type, celltype)
  p <- plot_cell_types(cell_type_banksy_counts, "celltype")
  ggsave(
  filename =file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_cluster_plot.png")),
  plot = p,
  width = 10,
  height = 7,
  units = "in")
  counts <- table(left0$first_type, left0@meta.data[[group_col]])
  
  # Normalize by column (so each group sums to 100%)
  props <- prop.table(counts, margin = 2) * 100
  
  # Define filename
  heatmap_file_png <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_heatmap.png"))
  heatmap_file_pdf <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_heatmap.pdf"))

  reds_palette <- colorRampPalette(c("white", "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15"))(100)
  

  annotation_col <- data.frame(
    Group = colnames(props)
  )
  rownames(annotation_col) <- colnames(props)

  names(annotation_col) <- group_col

  annotation_row <- data.frame(
    RCTD_Type = rownames(props)
  )
  rownames(annotation_row) <- rownames(props)

  get_colors <- function(n) {
    if(n <= 9) return(RColorBrewer::brewer.pal(max(3, n), "Set1")[1:n])
    if(n <= 12) return(RColorBrewer::brewer.pal(n, "Set3"))
    return(rainbow(n))
  }

  anno_colors <- list()
  anno_colors[[group_col]] <- setNames(get_colors(ncol(props)), colnames(props))
  anno_colors[["RCTD_Type"]] <- setNames(get_colors(nrow(props)), rownames(props))

  tryCatch({
    # Save as PDF
    pdf(file = heatmap_file_pdf, width = 10, height = 8)
    pheatmap(props, 
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
             main = paste0("RCTD first_type proportion by ", group_col))
    dev.off()
             
    message(paste0("Heatmap saved to ", heatmap_file_png))
  }, error = function(e) {
    message("Error generating heatmap: ", e$message)
    try({dev.off()}, silent = TRUE) # Ensure device is closed on error
  })
           
} else {
  warning(paste0("Column '", group_col, "' not found in spatial metadata. Skipping heatmap."))
  message("Available columns: ", paste(colnames(left0@meta.data), collapse = ", "))
}


if (all(c("spot_class", "first_type") %in% colnames(left0@meta.data))) {  
  # Prepare data: count frequency of each spot_class within each first_type
  spot_class_counts <- left0@meta.data %>%
    dplyr::select(first_type, spot_class) %>%
    dplyr::group_by(first_type, spot_class) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    dplyr::group_by(first_type) %>%
    dplyr::mutate(proportion = count / sum(count) * 100)
  
  # Define colors for spot_class (singlet, doublet_certain, doublet_uncertain, reject)
  # Using colors similar to the example plot
  spot_class_colors <- c(
    "singlet" = "#377EB8",           # Blue
    "doublet_certain" = "#FF7F00",   # Orange
    "doublet_uncertain" = "#984EA3", # Purple
    "reject" = "#E41A1C"             # Red
  )
  
  # Plot
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
  
  # Save plot
  bar_plot_file <- file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_spot_class_bar.png"))
  ggsave(filename = bar_plot_file, plot = p_bar, width = 12, height = 6, dpi = 300)
  message(paste0("Spot class bar plot saved to ", bar_plot_file))
}
message("RCTD Analysis Completed.")