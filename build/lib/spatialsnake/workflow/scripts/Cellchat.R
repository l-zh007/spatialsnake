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
  make_option(c("--interaction_length"), type="integer", default=150),
  make_option(c("--run_type"), type="character")
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
  vals <- trim_vec(unlist(strsplit(x, ",")))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) {
    return(NULL)
  }
  vals
}

normalize_run_type <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") {
    return("")
  }
  gsub("[[:space:]_-]+", "", tolower(as.character(x)))
}

is_visium_family <- function(run_type) {
  normalize_run_type(run_type) %in% c("visium", "visiumhd", "visiumsegment")
}

is_stereoseq_type <- function(run_type) {
  normalize_run_type(run_type) == "stereoseq"
}

is_merfish_type <- function(run_type) {
  normalize_run_type(run_type) %in% c("merfish", "merscope")
}

is_xenium_type <- function(run_type) {
  normalize_run_type(run_type) == "xenium"
}

infer_visium_hd_bin_size <- function(input_path) {
  patterns <- c("square_([0-9]+)um", "_([0-9]+)um")
  for (pattern in patterns) {
    hit <- regexpr(pattern, input_path, perl = TRUE)
    if (hit[1] > 0) {
      matched <- regmatches(input_path, hit)
      value <- sub(pattern, "\\1", matched, perl = TRUE)
      numeric_value <- suppressWarnings(as.numeric(value))
      if (!is.na(numeric_value) && numeric_value > 0) {
        return(numeric_value)
      }
    }
  }
  NA_real_
}

resolve_default_spot_size <- function(run_type, configured_spot_size, input_path) {
  rt <- normalize_run_type(run_type)
  if (rt == "visium") {
    return(65)
  }
  if (rt == "visiumhd") {
    inferred <- infer_visium_hd_bin_size(input_path)
    if (!is.na(inferred)) {
      return(inferred)
    }
    return(configured_spot_size)
  }
  if (rt %in% c("visiumsegment", "merfish", "merscope", "xenium")) {
    if (!is.null(configured_spot_size) && !is.na(configured_spot_size) && configured_spot_size != 65) {
      return(configured_spot_size)
    }
    return(10)
  }
  configured_spot_size
}

build_scale_factors <- function(spot_size, ratio) {
  if (is.null(ratio) || is.na(ratio) || ratio <= 0) {
    ratio <- 1
  }
  list(spot.diameter = spot_size, spot = spot_size / ratio)
}

build_spatial_factors <- function(ratios, spot_sizes, sample_names = NULL) {
  if (length(spot_sizes) == 1 && length(ratios) > 1) {
    spot_sizes <- rep(spot_sizes, length(ratios))
  }
  sf_df <- data.frame(ratio = ratios, tol = spot_sizes / 2)
  if (!is.null(sample_names) && length(sample_names) == nrow(sf_df)) {
    rownames(sf_df) <- sample_names
  }
  sf_df
}

parse_visium_json_spec <- function(spec, spot_size, run_type) {
  if (is.null(spec) || spec == "") {
    stop(paste0(run_type, " requires scalefactors_json.json in sample.txt third column"))
  }
  if (!file.exists(spec)) {
    stop(paste("scale_factors file not found:", spec))
  }
  print(spec)
  sf <- jsonlite::fromJSON(txt = spec)
  if (is.null(sf$spot_diameter_fullres) || is.na(sf$spot_diameter_fullres) || sf$spot_diameter_fullres <= 0) {
    stop(paste("spot_diameter_fullres missing or invalid in", spec))
  }
  ratio <- spot_size / sf$spot_diameter_fullres
  list(
    ratio = ratio,
    spot_size = spot_size,
    scale.factors = list(spot.diameter = spot_size, spot = sf$spot_diameter_fullres),
    source = spec
  )
}

parse_stereoseq_scale_spec <- function(spec, configured_spot_size) {
  if (is.null(spec) || spec == "") {
    stop("Stereo-seq requires sample.txt third column to be a bin_size number or 'cellbin'")
  }
  spec_norm <- normalize_run_type(spec)
  if (spec_norm %in% c("cellbin", "cell", "adjustedcellbin", "adjusted")) {
    spot_size <- if (!is.null(configured_spot_size) && !is.na(configured_spot_size) && configured_spot_size != 65) configured_spot_size else 10
    return(list(
      ratio = 0.5,
      spot_size = spot_size,
      scale.factors = build_scale_factors(spot_size, 0.5),
      source = spec
    ))
  }
  bin_size <- suppressWarnings(as.numeric(spec))
  if (is.na(bin_size) || bin_size <= 0) {
    stop(paste("Invalid Stereo-seq spec:", spec, "Expected a positive bin_size or 'cellbin'"))
  }
  spot_size <- bin_size * 0.5
  list(
    ratio = 0.5,
    spot_size = spot_size,
    scale.factors = build_scale_factors(spot_size, 0.5),
    source = spec
  )
}

resolve_platform_spatial_params <- function(run_type, input_path, scale_factor_spec, scale_factor_spec_list, sample_names, configured_spot_size) {
  rt <- normalize_run_type(run_type)
  default_spot_size <- resolve_default_spot_size(run_type, configured_spot_size, input_path)
  has_multi_specs <- !is.null(scale_factor_spec_list) && length(scale_factor_spec_list) > 0

  if (rt %in% c("visium", "visiumhd", "visiumsegment")) {
    if (has_multi_specs) {
      parsed <- lapply(scale_factor_spec_list, function(spec) parse_visium_json_spec(spec, default_spot_size, run_type))
      return(list(
        scale.factors = parsed[[1]]$scale.factors,
        conversion.factor = 1,
        spatial.factors = build_spatial_factors(
          vapply(parsed, function(item) item$ratio, numeric(1)),
          default_spot_size,
          sample_names
        ),
        spot_size = default_spot_size
      ))
    }
    parsed <- parse_visium_json_spec(scale_factor_spec, default_spot_size, run_type)
    return(list(
      scale.factors = parsed$scale.factors,
      conversion.factor = parsed$ratio,
      spatial.factors = build_spatial_factors(parsed$ratio, parsed$spot_size),
      spot_size = parsed$spot_size
    ))
  }

  if (rt == "stereoseq") {
    if (has_multi_specs) {
      parsed <- lapply(scale_factor_spec_list, function(spec) parse_stereoseq_scale_spec(spec, configured_spot_size))
      return(list(
        scale.factors = parsed[[1]]$scale.factors,
        conversion.factor = 1,
        spatial.factors = build_spatial_factors(
          vapply(parsed, function(item) item$ratio, numeric(1)),
          vapply(parsed, function(item) item$spot_size, numeric(1)),
          sample_names
        ),
        spot_size = parsed[[1]]$spot_size
      ))
    }
    parsed <- parse_stereoseq_scale_spec(scale_factor_spec, configured_spot_size)
    return(list(
      scale.factors = parsed$scale.factors,
      conversion.factor = parsed$ratio,
      spatial.factors = build_spatial_factors(parsed$ratio, parsed$spot_size),
      spot_size = parsed$spot_size
    ))
  }

  if (rt %in% c("merfish", "merscope", "xenium")) {
    spot_size <- resolve_default_spot_size(run_type, configured_spot_size, input_path)
    ratios <- if (!is.null(sample_names) && length(sample_names) > 0) rep(1, length(sample_names)) else 1
    return(list(
      scale.factors = build_scale_factors(spot_size, 1),
      conversion.factor = 1,
      spatial.factors = build_spatial_factors(ratios, spot_size, sample_names),
      spot_size = spot_size
    ))
  }

  scale.factors <- build_scale_factors(configured_spot_size, 1)
  list(
    scale.factors = scale.factors,
    conversion.factor = 1,
    spatial.factors = build_spatial_factors(1, configured_spot_size),
    spot_size = configured_spot_size
  )
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


data <- tryCatch({
  Seurat::GetAssayData(obj, layer = "data")
}, error = function(e_data) {
  tryCatch({
    obj <- JoinLayers(obj)
    Seurat::GetAssayData(obj, layer = "data")
  }, error = function(e_joined_data) {
    tryCatch({
      Seurat::GetAssayData(obj, layer = "count")
    }, error = function(e_count) {
      obj <- JoinLayers(obj)
      Seurat::GetAssayData(obj, layer = "count")
    })
  })
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
  scale_factors_list <- split_vec(opt$scale_factors_list)
  if (!is.null(scale_factors_list) && !is.null(sample_names) && length(scale_factors_list) != length(sample_names)) {
    stop("scale_factors_list and sample_names length mismatch")
  }

  platform_params <- resolve_platform_spatial_params(
    run_type = opt$run_type,
    input_path = opt$input_path,
    scale_factor_spec = opt$scale_factors,
    scale_factor_spec_list = scale_factors_list,
    sample_names = sample_names,
    configured_spot_size = opt$spot_size
  )
  scale.factors <- platform_params$scale.factors
  conversion.factor <- platform_params$conversion.factor
  spatial.factors <- platform_params$spatial.factors
  message(paste("Resolved run_type:", opt$run_type))
  message(paste("Resolved spot_size:", platform_params$spot_size))
  message(paste("Conversion factor:", conversion.factor))
  cellchat <- createCellChat(object = data, meta = meta, group.by = "celltypes", datatype = "spatial", coordinates = spatial.mat, scale.factors = scale.factors)
  
  # Ensure scale.distance satisfies CellChat's requirement (min scaled distance >= 1)
  if (nrow(spatial.mat) > 1) {
    coords_sorted_x <- spatial.mat[order(spatial.mat[,1], spatial.mat[,2]), , drop = FALSE]
    dx <- diff(coords_sorted_x[,1])
    dy <- diff(coords_sorted_x[,2])
    dists_x <- sqrt(dx^2 + dy^2)
    
    coords_sorted_y <- spatial.mat[order(spatial.mat[,2], spatial.mat[,1]), , drop = FALSE]
    dx2 <- diff(coords_sorted_y[,1])
    dy2 <- diff(coords_sorted_y[,2])
    dists_y <- sqrt(dx2^2 + dy2^2)
    
    dists_valid <- c(dists_x[dists_x > 0], dists_y[dists_y > 0])
    if (length(dists_valid) > 0) {
      min_dist <- min(dists_valid)
      scaled_dist <- conversion.factor * min_dist
      if (scaled_dist < 1 || scaled_dist > 2) {
        conversion.factor <- 1.5 / min_dist
        message(paste("Adjusted conversion.factor to", conversion.factor, "to satisfy CellChat distance requirements"))
      }
    }
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
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_count_heatmap.png")), width = 1200, height = 900)
p_count <- netVisual_heatmap(cellchat, measure = "count")
if (!is.null(p_count)) {
  print(p_count)
}
dev.off()
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_heatmap.png")), width = 1200, height = 900)
p_weight <- netVisual_heatmap(cellchat, measure = "weight")
if (!is.null(p_weight)) {
  print(p_weight)
}
dev.off()
pathway_auto <- NULL
if (!is.null(cellchat@netP$pathways) && length(cellchat@netP$pathways) > 0) {
  pathway_auto <- cellchat@netP$pathways[1]
}
cellchat_role <- tryCatch({
  netAnalysis_computeCentrality(cellchat, slot.name = "netP")
}, error = function(e) {
  NULL
})
if (!is.null(cellchat_role) && !is.null(pathway_auto)) {
  png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_signaling_role_network.png")), width = 1500, height = 700)
  p_role_network <- tryCatch({
    netAnalysis_signalingRole_network(cellchat_role, width = 15, height = 6, font.size = 10)
  }, error = function(e) {
    NULL
  })
  if (!is.null(p_role_network)) {
    print(p_role_network)
  }
  dev.off()
  png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_signaling_role_scatter.png")), width = 1000, height = 800)
  p_role_scatter <- tryCatch({
    netAnalysis_signalingRole_scatter(cellchat_role)
  }, error = function(e) {
    NULL
  })
  if (!is.null(p_role_scatter)) {
    print(p_role_scatter)
  }
  dev.off()
  png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_signaling_role_outgoing.png")), width = 1200, height = 900)
  p_role_out <- tryCatch({
    netAnalysis_signalingRole_heatmap(cellchat_role, pattern = "outgoing")
  }, error = function(e) {
    NULL
  })
  if (!is.null(p_role_out)) {
    print(p_role_out)
  }
  dev.off()
  png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_signaling_role_incoming.png")), width = 1200, height = 900)
  p_role_in <- tryCatch({
    netAnalysis_signalingRole_heatmap(cellchat_role, pattern = "incoming")
  }, error = function(e) {
    NULL
  })
  if (!is.null(p_role_in)) {
    print(p_role_in)
  }
  dev.off()
}
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
