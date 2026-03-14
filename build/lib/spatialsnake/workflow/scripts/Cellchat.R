library(optparse)
library(Seurat)
library(CellChat)
library(jsonlite)
library(Matrix)
library(dplyr)
library(ggplot2)
library(future)
library(schard)
option_list <- list(
  make_option(c("-s", "--input_path"), type="character"),
  make_option(c("-o", "--output_dir"), type="character"),
  make_option(c("--sample_id"), type="character", default="sample"),
  make_option(c("--celltype_col"), type="character", default="celltype"),
  make_option(c("--assay"), type="character", default="Spatial"),
  make_option(c("--species"), type="character", default="human"),
  make_option(c("--min_cells"), type="integer", default=10),
  make_option(c("--nworkers"), type="integer", default=4),
  make_option(c("--is_single_cell"), type="logical"),
  make_option(c("--scale_factors"), type="character", default=""),
  make_option(c("--scale_factors_list"), type="character", default=""),
  make_option(c("--sample_names"), type="character", default=""),
  make_option(c("--spot_size"), type="double", default=65),
  make_option(c("--trim"), type="double", default=0.1),
  make_option(c("--interaction_length"), type="integer", default=150)
)
opt <- parse_args(OptionParser(option_list=option_list))
print(opt$is_single_cell)
print(opt$scale_factors_list)
print(opt$sample_names)
print(opt$input_path)

trim_vec <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

split_vec <- function(x) {
  if (is.null(x) || x == "" || tolower(x) %in% c("null", "none")) {
    return(NULL)
  }
  trim_vec(unlist(strsplit(x, ",")))
}

h5ad_spatial_loading <- function(input_h5ad) {
  obj <- tryCatch({
    obj <- schard::h5ad2seurat_spatial(input_h5ad, use.raw = TRUE)
    obj <- NormalizeData(obj)
  }, error = function(e) {
    message(paste("Error normalizing spatial data:", e$message))
    obj <- schard::h5ad2seurat_spatial(input_h5ad, use.raw = FALSE)
  })
  return(obj)
}
load_seurat <- function(file_path, is_single_cell = opt$is_single_cell) {
  if (!is_single_cell) {
    if (grepl("\\.h5ad$", file_path, ignore.case = TRUE)) {
      message(paste("Loading spatial data from", file_path))
      h5ad_spatial_loading(file_path)
    } else {
      readRDS(file_path)
    }
  } else {
    if (grepl("\\.rds$", file_path, ignore.case = TRUE)) {
      readRDS(file_path)
    } else if (grepl("\\.h5ad$", file_path, ignore.case = TRUE)) {
      tryCatch({
        schard::h5ad2seurat(file_path, use.raw = TRUE)
      }, error = function(e) {
        schard::h5ad2seurat(file_path, use.raw = FALSE)
      })
    } else {
      stop("unsupported file")
    }
  }
}


obj <- load_seurat(opt$input_path, is_single_cell = opt$is_single_cell)


data <-tryCatch({Seurat::GetAssayData(obj, layer = "count")
}, error = function(e) {
        obj <- JoinLayers(obj)
        Seurat::GetAssayData(obj, layer = "count")
})


Idents(obj) <- obj@meta.data[[opt$celltype_col]]
meta <- data.frame(celltypes = Idents(obj), row.names = names(Idents(obj)))
CellChatDB <- if (tolower(opt$species) == "mouse") CellChatDB.mouse else CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")




if (opt$is_single_cell) {
  cellchat <- createCellChat(object = data, meta = meta, group.by = "celltypes", datatype = "scRNAseq")
} else {
  img_names <- Images(obj)
  sample_names <- split_vec(opt$sample_names)
  img_name <- if (length(img_names) > 0) img_names[1] else NULL
  if ("boundaries" %in% slotNames(obj) && !is.null(tryCatch(obj@boundaries$centroids, error=function(e) NULL))) {
    spatial.locs <- obj@boundaries$centroids@coords
  } else {
    if (length(img_names) > 1) {
      img_order <- if (!is.null(sample_names)) intersect(sample_names, img_names) else img_names
      if (length(img_order) == 0) {
        img_order <- img_names
      }
      spatial_list <- lapply(img_order, function(img) {
        Seurat::GetTissueCoordinates(obj, image = img, cols = c("imagerow", "imagecol"))
      })
      spatial.locs <- do.call(rbind, spatial_list)
    } else {
      spatial.locs <- Seurat::GetTissueCoordinates(obj, image = img_name,cols = c("imagerow", "imagecol"))
    }
  }
  spatial.df <- data.frame(x = spatial.locs$imagerow, y = spatial.locs$imagecol)
  rownames(spatial.df) <- rownames(spatial.locs)
  common_cells <- intersect(rownames(spatial.df), colnames(data))
  spatial.df <- spatial.df[common_cells, , drop = FALSE]
  data <- data[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  spatial.mat <- as.matrix(spatial.df)
  scale.factors <- NULL
  conversion.factor <- 1
  spatial.factors <- NULL
  if (opt$scale_factors != "" && file.exists(opt$scale_factors)) {
    sf <- jsonlite::fromJSON(txt = opt$scale_factors)
    if (!is.null(sf$spot_diameter_fullres)) {
      spot_diameter <- sf$spot_diameter_fullres
      scale.factors <- list(spot.diameter = spot_diameter, spot = spot_diameter)
      conversion.factor <- spot_diameter / sf$spot_diameter_fullres
    } else if (!is.null(sf$spot_diameter)) {
      spot_diameter <- sf$spot_diameter
      scale.factors <- list(spot.diameter = spot_diameter, spot = spot_diameter)
      conversion.factor <- 1
    }
  }

  scale_factors_list <- split_vec(opt$scale_factors_list)
  if (!is.null(scale_factors_list) || !is.null(sample_names)) {
    if (is.null(scale_factors_list) || is.null(sample_names)) {
      stop("scale_factors_list and sample_names must be provided together")
    }
    if (length(scale_factors_list) != length(sample_names)) {
      stop("scale_factors_list and sample_names length mismatch")
    }
    ratios <- c()
    for (i in seq_along(scale_factors_list)) {
      sf_file <- scale_factors_list[[i]]
      if (!file.exists(sf_file)) {
        stop(paste("scale_factors file not found:", sf_file))
      }
      sf_item <- jsonlite::fromJSON(txt = sf_file)
      if (is.null(sf_item$spot_diameter_fullres)) {
        stop(paste("spot_diameter_fullres missing in", sf_file))
      }
      ratios <- c(ratios, opt$spot_size / sf_item$spot_diameter_fullres)
    }
    spatial.factors <- data.frame(ratio = ratios, tol = opt$spot_size / 2)
    rownames(spatial.factors) <- sample_names
    conversion.factor <- 1
  }

  message(paste("Conversion factor:", conversion.factor))
  if (is.null(scale.factors)) {
    cellchat <- createCellChat(object = data, meta = meta, group.by = "celltypes", datatype = "spatial", coordinates = spatial.mat)
  } else {
    cellchat <- createCellChat(object = data, meta = meta, group.by = "celltypes", datatype = "spatial", coordinates = spatial.mat, scale.factors = scale.factors)
  }
}
cellchat@DB <- CellChatDB.use
options(future.globals.maxSize = 8 * 1024^3)
future::plan("multicore", workers = opt$nworkers)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
if (opt$is_single_cell) {
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = opt$trim)
} else {
  use_spatial_factors <- !is.null(spatial.factors)
  if (use_spatial_factors && "spatial.factors" %in% names(formals(CellChat::computeCommunProb))) {
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = opt$trim, distance.use = TRUE, scale.distance = conversion.factor, interaction.length = opt$interaction_length, spatial.factors = spatial.factors)
  } else {
    if (use_spatial_factors) {
      message("spatial.factors is ignored because computeCommunProb does not accept it")
    }
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = opt$trim, distance.use = TRUE, scale.distance = conversion.factor, interaction.length = opt$interaction_length)
  }
}
cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
cellchat <- computeCommunProbPathway(cellchat)
net <- subsetCommunication(cellchat)
netP <- subsetCommunication(cellchat, slot.name = "netP")
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(opt$output_dir, "cellchat.rds"))
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_network.png")), width = 1600, height = 800)
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()
pdf(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_network.pdf")), width = 14, height = 7)
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_infoflow_bar.png")), width = 1200, height = 900)
rankNet(cellchat, mode = "single", stacked = TRUE, do.stat = FALSE)
dev.off()
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_heatmap.png")), width = 1200, height = 900)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()
counts <- cellchat@net$count
weights <- cellchat@net$weight
df_counts <- as.data.frame(counts)
df_counts$source <- rownames(df_counts)
df_counts <- tidyr::pivot_longer(df_counts, cols = -source, names_to = "target", values_to = "count")
df_weights <- as.data.frame(weights)
df_weights$source <- rownames(df_weights)
df_weights <- tidyr::pivot_longer(df_weights, cols = -source, names_to = "target", values_to = "weight")
stats <- dplyr::left_join(df_counts, df_weights, by = c("source","target"))
net_lr <- net
net_lr$sample_id <- opt$sample_id
if (all(c("ligand","receptor") %in% names(net_lr))) {
  net_lr$lr_pair <- paste(net_lr$ligand, net_lr$receptor, sep = "_")
}
pathway_summary <- NULL
pathway_pairs <- NULL
if ("pathway_name" %in% names(net_lr)) {
  pathway_summary <- net_lr %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(
      interaction_n = dplyr::n(),
      lr_n = dplyr::n_distinct(dplyr::coalesce(interaction_name, interaction_name_2, lr_pair)),
      source_n = dplyr::n_distinct(source),
      target_n = dplyr::n_distinct(target),
      prob_sum = sum(prob, na.rm = TRUE),
      prob_mean = mean(prob, na.rm = TRUE),
      prob_max = max(prob, na.rm = TRUE),
      pval_min = if ("pval" %in% names(net_lr)) min(pval, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  pathway_pairs <- net_lr %>%
    dplyr::group_by(pathway_name, source, target) %>%
    dplyr::summarise(
      interaction_n = dplyr::n(),
      prob_sum = sum(prob, na.rm = TRUE),
      prob_mean = mean(prob, na.rm = TRUE),
      pval_min = if ("pval" %in% names(net_lr)) min(pval, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}
lr_summary <- NULL
if ("lr_pair" %in% names(net_lr)) {
  lr_summary <- net_lr %>%
    dplyr::group_by(lr_pair) %>%
    dplyr::summarise(
      interaction_n = dplyr::n(),
      pathway_n = if ("pathway_name" %in% names(net_lr)) dplyr::n_distinct(pathway_name) else NA_integer_,
      prob_sum = sum(prob, na.rm = TRUE),
      prob_mean = mean(prob, na.rm = TRUE),
      prob_max = max(prob, na.rm = TRUE),
      pval_min = if ("pval" %in% names(net_lr)) min(pval, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}
write.csv(stats, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_stats.csv")), row.names = FALSE)
write.csv(net_lr, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_lr.csv")), row.names = FALSE)
if (!is.null(lr_summary)) {
  write.csv(lr_summary, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_lr_summary.csv")), row.names = FALSE)
}
if (!is.null(pathway_pairs)) {
  write.csv(pathway_pairs, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_pathway_pairs.csv")), row.names = FALSE)
}
if (!is.null(pathway_summary)) {
  write.csv(pathway_summary, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_pathway_summary.csv")), row.names = FALSE)
}
if (is.data.frame(netP) || is.matrix(netP)) {
  write.csv(netP, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_pathway_net.csv")), row.names = FALSE)
}
