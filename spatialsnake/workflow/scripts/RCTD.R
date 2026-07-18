library(optparse)
library(Seurat)
library(spacexr)
library(Matrix)

option_list <- list(
  make_option(c("-s", "--spatial_input"), type="character", help="Path to spatial data (Seurat .rds or .h5Seurat or .h5ad)"),
  make_option(c("-r", "--sc_input"), type="character", help="Path to single cell data (Seurat .rds or .h5Seurat or .h5ad)"),
  make_option(c("-o", "--output_dir"), type="character", help="Directory to save results"),
  make_option(c("--sample_id"), type="character", default="sample", help="Sample ID"),
  make_option(c("--mode"), type="character", default="doublet", help="RCTD mode: 'doublet' or 'full'"),
  make_option(c("--max_cores"), type="integer", default=8, help="Number of cores for RCTD"),
  make_option(c("--sc_cell_type_col"), type="character", default="celltype", help="Column name for cell types in SC data"),
  make_option(c("--spatial_cell_type_col"), type="character", default="celltype", help="Column name for cell types or clusters in spatial metadata"),
  make_option(c("--group_by"), type="character", default="sample", help="Column in Spatial Seurat object to group by (X-axis) for heatmap")
)

opt <- parse_args(OptionParser(option_list=option_list))
opt$mode <- tolower(trimws(opt$mode))
if (!opt$mode %in% c("doublet", "full")) {
  stop("--mode must be either 'doublet' or 'full'.")
}
if (is.na(opt$max_cores) || opt$max_cores < 1L) {
  stop("--max_cores must be a positive integer.")
}
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
    sc_counts <- GetAssayData(merge_seurat_object, assay = "RNA", layer = "counts")
  }
} else {
  # 尝试使用默认 assay
  sc_counts <- GetAssayData(merge_seurat_object, layer = "counts")
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
  st_count <- GetAssayData(left0, layer = "counts")
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

if (!identical(ncol(st_loc), 2L)) {
  if (all(c("x", "y") %in% colnames(st_loc))) {
    st_loc <- st_loc[, c("x", "y"), drop = FALSE]
  } else if (ncol(st_loc) >= 2L) {
    st_loc <- st_loc[, seq_len(2), drop = FALSE]
  } else {
    stop("Spatial coordinates must have at least two columns after loading.")
  }
}

colnames(st_loc) <- c("x", "y")

if (ncol(st_loc) != 2L) {
  stop("Spatial coordinates normalization failed: coords must contain exactly two columns named x and y.")
}

if (length(colnames(left0)) == nrow(st_loc)) {
  rownames(st_loc) <- colnames(left0)
}

message(
  "Coordinate table after normalization: ",
  nrow(st_loc), " x ", ncol(st_loc),
  " [", paste(colnames(st_loc), collapse = ", "), "]"
)

# 创建 SpatialRNA 对象
query <- spacexr::SpatialRNA(st_loc, st_count)

# 5. 运行 RCTD
message("Running RCTD (mode: ", opt$mode, ")...")
RCTD <- create.RCTD(query, reference, max_cores = opt$max_cores)
RCTD <- run.RCTD(RCTD, doublet_mode = opt$mode)

# 6. 保存结果。绘图统一交由 Python 后处理，避免 R/Python 重复实现。
message("Saving results to ", opt$output_dir)

raw_weights <- RCTD@results$weights
if (is.null(raw_weights) || length(dim(raw_weights)) != 2L || nrow(raw_weights) == 0L || ncol(raw_weights) == 0L) {
  stop("RCTD did not return a non-empty weights matrix.")
}

# full 模式的原始权重不保证逐行和为 1；官方建议使用 normalize_weights。
norm_weights <- as.matrix(normalize_weights(raw_weights))
storage.mode(norm_weights) <- "double"

if (is.null(rownames(norm_weights)) || anyNA(rownames(norm_weights)) || anyDuplicated(rownames(norm_weights))) {
  stop("RCTD weights must have unique, non-missing spot identifiers.")
}
if (is.null(colnames(norm_weights)) || anyNA(colnames(norm_weights)) ||
    any(trimws(colnames(norm_weights)) == "") || anyDuplicated(colnames(norm_weights))) {
  stop("RCTD weights must have unique, non-empty cell-type names.")
}
if (any(!is.finite(norm_weights))) {
  stop("Normalized RCTD weights contain non-finite values.")
}
if (any(norm_weights < -1e-8)) {
  stop("Normalized RCTD weights contain negative values.")
}

# 容忍浮点误差，并再次保证导出的每个有效 spot 权重和精确为 1。
norm_weights[norm_weights < 0] <- 0
weight_sums <- rowSums(norm_weights)
if (any(!is.finite(weight_sums)) || any(weight_sums <= 0)) {
  stop("Normalized RCTD weights contain rows with non-positive totals.")
}
norm_weights <- norm_weights / weight_sums

results_df <- RCTD@results$results_df
if (is.null(results_df)) {
  # full 模式官方仅返回 weights；结果 CSV 只保留拟合后的 spot ID。
  write.csv(
    data.frame(spot_id = rownames(norm_weights), stringsAsFactors = FALSE),
    file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_results.csv")),
    row.names = FALSE
  )
} else {
  results_df <- as.data.frame(results_df)
  if (is.null(rownames(results_df)) || anyDuplicated(rownames(results_df))) {
    stop("RCTD results_df must have unique spot identifiers.")
  }
  write.csv(
    results_df,
    file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_results.csv")),
    row.names = TRUE
  )
}

write.csv(
  norm_weights,
  file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_weights.csv")),
  row.names = TRUE
)

if (!is.null(results_df)) {
  left0 <- AddMetaData(left0, metadata = results_df)
}

saveRDS(left0, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD_seurat.rds")))
saveRDS(RCTD, file = file.path(opt$output_dir, paste0(opt$sample_id, "_RCTD.rds")))

message("RCTD analysis completed. Visualization will be generated from the exported tables in Python.")
