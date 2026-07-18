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
  make_option(c("--db_subset"), type="character", default="all_interactions"),
  make_option(c("--trim"), type="double", default=0.1),
  make_option(c("--interaction_length"), type="integer", default=250),
  make_option(c("--run_type"), type="character"),
  make_option(c("--pathways"), type="character", default=""),
  make_option(c("--top_pathways"), type="integer", default=3),
  make_option(c("--focus_cells"), type="character", default=""),
  make_option(c("--source_cells"), type="character", default=""),
  make_option(c("--target_cells"), type="character", default=""),
  make_option(c("--pair_lr_use"), type="character", default=""),
  make_option(c("--cell_pairs"), type="character", default=""),
  make_option(c("--lr_pairs"), type="character", default=""),
  make_option(c("--top_cell_pairs"), type="integer", default=3),
  make_option(c("--bubble_top_lr"), type="integer", default=20),
  make_option(c("--plot_advanced"), type="logical", default=TRUE),
  make_option(c("--future_max_size_gb"), type="double", default=64)
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
message(paste("CellChat input:", opt$input_path))
message(paste("CellChat run_type:", opt$run_type))
message(paste("CellChat mode:", ifelse(isTRUE(opt$is_single_cell), "single-cell", "spatial")))
message(paste("CellChat interaction_length:", opt$interaction_length))
message(paste("CellChat configured spot_size:", opt$spot_size))
message(paste("CellChat scale_factors_list:", opt$scale_factors_list))
message(paste("CellChat sample_names:", opt$sample_names))
message(paste("CellChat future.globals.maxSize (GB):", opt$future_max_size_gb))

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

valid_named_values <- function(requested, available, label) {
  values <- split_vec(requested)
  if (is.null(values)) {
    return(NULL)
  }
  matched <- intersect(values, available)
  missing <- setdiff(values, available)
  if (length(missing) > 0) {
    warning(paste("Requested", label, "not found and will be ignored:", paste(missing, collapse = ", ")))
  }
  if (length(matched) == 0) {
    return(NULL)
  }
  matched
}

select_pathways <- function(requested, available, pathway_summary = NULL, top_n = 5) {
  if (is.null(available) || length(available) == 0) {
    return(character(0))
  }
  requested_pathways <- valid_named_values(requested, available, "CellChat pathway")
  if (!is.null(requested_pathways) && length(requested_pathways) > 0) {
    return(requested_pathways)
  }
  top_n <- max(1, top_n)
  if (!is.null(pathway_summary) && "pathway_name" %in% names(pathway_summary) && "prob_sum" %in% names(pathway_summary)) {
    ordered <- pathway_summary %>%
      dplyr::arrange(dplyr::desc(prob_sum)) %>%
      dplyr::pull(pathway_name)
    ordered <- ordered[ordered %in% available]
    if (length(ordered) > 0) {
      return(head(unique(ordered), top_n))
    }
  }
  head(available, top_n)
}

normalise_lr_pair <- function(x) {
  gsub("[|: /-]+", "_", trim_vec(x))
}

split_pair_token <- function(x) {
  token <- trim_vec(x)
  if (grepl("<->", token, fixed = TRUE)) {
    parts <- trim_vec(unlist(strsplit(token, "<->", fixed = TRUE)))
    if (length(parts) == 2 && all(nzchar(parts))) {
      return(list(c(parts[1], parts[2]), c(parts[2], parts[1])))
    }
  }
  if (grepl("->", token, fixed = TRUE)) {
    parts <- trim_vec(unlist(strsplit(token, "->", fixed = TRUE)))
  } else if (grepl("\\|", token)) {
    parts <- trim_vec(unlist(strsplit(token, "\\|")))
  } else {
    parts <- character(0)
  }
  if (length(parts) == 2 && all(nzchar(parts))) {
    return(list(c(parts[1], parts[2])))
  }
  list()
}

parse_cell_pairs <- function(x) {
  values <- split_vec(x)
  if (is.null(values)) {
    return(NULL)
  }
  pairs <- unlist(lapply(values, split_pair_token), recursive = FALSE)
  if (length(pairs) == 0) {
    return(NULL)
  }
  unique(vapply(pairs, function(pair) paste(pair[1], pair[2], sep = "|"), character(1)))
}

parse_lr_pairs <- function(x) {
  values <- split_vec(x)
  if (is.null(values)) {
    return(NULL)
  }
  unique(normalise_lr_pair(values))
}

filter_lr_pairs <- function(df, requested_pairs) {
  if (is.null(requested_pairs) || nrow(df) == 0) {
    return(df)
  }
  pair_match <- rep(FALSE, nrow(df))
  if ("lr_pair" %in% names(df)) {
    pair_match <- pair_match | normalise_lr_pair(df$lr_pair) %in% requested_pairs
  }
  if ("interaction_name" %in% names(df)) {
    pair_match <- pair_match | normalise_lr_pair(df$interaction_name) %in% requested_pairs
  }
  if ("interaction_name_2" %in% names(df)) {
    pair_match <- pair_match | normalise_lr_pair(df$interaction_name_2) %in% requested_pairs
  }
  df[pair_match, , drop = FALSE]
}

plot_dim <- function(n, base = 6, per = 0.28, min_size = 5, max_size = 16) {
  min(max_size, max(min_size, base + per * max(0, n)))
}

safe_label <- function(x) {
  label <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  gsub("^_+|_+$", "", label)
}

make_pair_lr_use <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !("interaction_name" %in% names(df))) {
    return(NULL)
  }
  data.frame(interaction_name = unique(as.character(df$interaction_name)))
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

is_default_rownames <- function(rn, n) {
  is.null(rn) || length(rn) != n || identical(rn, as.character(seq_len(n)))
}

pick_existing_column <- function(df, candidates) {
  if (is.null(df) || ncol(df) == 0) {
    return(NULL)
  }
  lower_names <- tolower(colnames(df))
  for (candidate in candidates) {
    hit <- which(lower_names == tolower(candidate))
    if (length(hit) > 0) {
      return(colnames(df)[hit[1]])
    }
  }
  NULL
}

pick_coordinate_columns <- function(df) {
  candidate_pairs <- list(
    c("imagerow", "imagecol"),
    c("x", "y"),
    c("row", "col"),
    c("pxl_row_in_fullres", "pxl_col_in_fullres"),
    c("array_row", "array_col"),
    c("x_centroid", "y_centroid"),
    c("centroid_x", "centroid_y"),
    c("spatial_1", "spatial_2"),
    c("spatial1", "spatial2")
  )
  lower_names <- tolower(colnames(df))
  for (pair in candidate_pairs) {
    idx <- match(tolower(pair), lower_names)
    if (all(!is.na(idx))) {
      return(colnames(df)[idx])
    }
  }
  NULL
}

standardize_coordinate_table <- function(coord_df, source_label = "spatial", data_cell_names = NULL) {
  coord_df <- as.data.frame(coord_df)
  if (nrow(coord_df) == 0) {
    stop(paste("No spatial coordinates were found for", source_label))
  }
  coord_cols <- pick_coordinate_columns(coord_df)
  if (is.null(coord_cols)) {
    stop(paste(
      "Unable to identify spatial coordinate columns for", source_label,
      "Available columns:", paste(colnames(coord_df), collapse = ", ")
    ))
  }

  rn <- rownames(coord_df)
  if (!is.null(data_cell_names) && length(intersect(rn, data_cell_names)) == 0) {
    cell_col <- pick_existing_column(coord_df, c("cell", "barcode", "barcodes", "cell_id", "spot_id"))
    if (!is.null(cell_col) && length(intersect(as.character(coord_df[[cell_col]]), data_cell_names)) > 0) {
      rn <- as.character(coord_df[[cell_col]])
    }
  } else if (is_default_rownames(rn, nrow(coord_df))) {
    cell_col <- pick_existing_column(coord_df, c("cell", "barcode", "barcodes", "cell_id", "spot_id"))
    if (!is.null(cell_col)) {
      rn <- as.character(coord_df[[cell_col]])
    }
  }

  if (is.null(rn) || length(rn) != nrow(coord_df) || any(is.na(rn)) || any(!nzchar(rn))) {
    stop(paste("Spatial coordinate row names are invalid for", source_label))
  }

  out <- data.frame(
    x = suppressWarnings(as.numeric(coord_df[[coord_cols[1]]])),
    y = suppressWarnings(as.numeric(coord_df[[coord_cols[2]]])),
    library_id = source_label,
    row.names = rn,
    stringsAsFactors = FALSE
  )
  finite <- is.finite(out$x) & is.finite(out$y)
  if (!all(finite)) {
    removed <- sum(!finite)
    warning(paste("Removed", removed, "spots/cells with non-finite spatial coordinates from", source_label))
    out <- out[finite, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    stop(paste("No finite spatial coordinates remain for", source_label))
  }
  out
}

extract_spatial_coordinates <- function(obj, sample_names = NULL, data_cell_names = NULL) {
  img_names <- Images(obj)
  if (length(img_names) > 0) {
    img_order <- img_names
    if (!is.null(sample_names) && length(sample_names) > 0) {
      matched <- intersect(sample_names, img_names)
      missing <- setdiff(sample_names, img_names)
      if (length(missing) > 0) {
        warning(paste("The following sample names were not found among Seurat images:", paste(missing, collapse = ", ")))
      }
      img_order <- c(matched, setdiff(img_names, matched))
    }
    spatial_list <- lapply(img_order, function(img) {
      coords <- Seurat::GetTissueCoordinates(obj, image = img)
      standardize_coordinate_table(coords, source_label = img, data_cell_names = data_cell_names)
    })
    spatial.df <- do.call(rbind, spatial_list)
  } else if ("boundaries" %in% slotNames(obj) && "centroids" %in% names(obj@boundaries) &&
             "coords" %in% slotNames(obj@boundaries$centroids)) {
    spatial.df <- standardize_coordinate_table(
      obj@boundaries$centroids@coords,
      source_label = "boundaries",
      data_cell_names = data_cell_names
    )
  } else if ("spatial" %in% Reductions(obj)) {
    emb <- Embeddings(obj, reduction = "spatial")
    if (ncol(emb) < 2) {
      stop("The spatial reduction must contain at least two coordinate dimensions")
    }
    spatial.df <- data.frame(
      x = as.numeric(emb[, 1]),
      y = as.numeric(emb[, 2]),
      library_id = "spatial",
      row.names = rownames(emb),
      stringsAsFactors = FALSE
    )
  } else {
    spatial.df <- standardize_coordinate_table(
      obj@meta.data,
      source_label = "metadata",
      data_cell_names = data_cell_names
    )
  }

  if (any(duplicated(rownames(spatial.df)))) {
    duplicated_cells <- unique(rownames(spatial.df)[duplicated(rownames(spatial.df))])
    stop(paste(
      "Spatial coordinates contain duplicated cell/spot identifiers. Examples:",
      paste(head(duplicated_cells, 5), collapse = ", ")
    ))
  }
  spatial.df
}

offset_coordinates_by_library <- function(spatial.df, interaction_length = 250) {
  libraries <- unique(as.character(spatial.df$library_id))
  libraries <- libraries[!is.na(libraries) & nzchar(libraries)]
  if (length(libraries) <= 1) {
    return(spatial.df)
  }
  x_range <- diff(range(spatial.df$x, na.rm = TRUE))
  y_range <- diff(range(spatial.df$y, na.rm = TRUE))
  gap <- max(c(x_range, y_range, interaction_length * 10, 1000), na.rm = TRUE) + interaction_length * 10
  for (i in seq_along(libraries)) {
    idx <- as.character(spatial.df$library_id) == libraries[i]
    spatial.df$x[idx] <- spatial.df$x[idx] + (i - 1) * gap
  }
  message(paste(
    "Applied library-wise coordinate offsets for", length(libraries),
    "spatial libraries because this CellChat version does not support spatial.factors."
  ))
  spatial.df
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

if (!(opt$celltype_col %in% colnames(obj@meta.data))) {
  stop(paste("celltype_col not found in Seurat metadata:", opt$celltype_col))
}
celltype_values <- as.character(obj@meta.data[[opt$celltype_col]])
celltype_values <- celltype_values[!is.na(celltype_values) & nzchar(celltype_values)]
celltype_counts <- table(celltype_values)
if (length(celltype_counts) < 2) {
  stop(paste("CellChat requires at least two valid cell groups in", opt$celltype_col))
}
small_groups <- names(celltype_counts)[celltype_counts < opt$min_cells]
if (length(small_groups) > 0) {
  warning(paste(
    "The following groups have fewer cells/spots than cellchat_min_cells and may be filtered:",
    paste(paste0(small_groups, "=", as.integer(celltype_counts[small_groups])), collapse = ", ")
  ))
}
message(paste("CellChat cell groups:", paste(names(celltype_counts), collapse = ", ")))
message(paste("CellChat min_cells:", opt$min_cells))

Idents(obj) <- obj@meta.data[[opt$celltype_col]]
meta <- data.frame(celltypes = Idents(obj), row.names = names(Idents(obj)))
CellChatDB <- if (tolower(opt$species) == "mouse") CellChatDB.mouse else CellChatDB.human
db_subset <- tolower(trimws(opt$db_subset))
CellChatDB.use <- switch(
  db_subset,
  "secreted_only" = subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation"),
  "secreted_ecm" = subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor"), key = "annotation"),
  "all_interactions" = CellChatDB,
  CellChatDB
)
message(paste("Using CellChat database subset:", db_subset))
message(paste("CellChat species:", opt$species))




if (opt$is_single_cell) {
  cellchat <- createCellChat(object = data, meta = meta, group.by = "celltypes", datatype = "RNA")
} else {
  sample_names <- split_vec(opt$sample_names)
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
  message("Resolved spatial.factors:")
  print(spatial.factors)

  spatial.df <- extract_spatial_coordinates(obj, sample_names = sample_names, data_cell_names = colnames(data))
  common_cells <- intersect(rownames(spatial.df), colnames(data))
  if (length(common_cells) == 0) {
    stop("No overlap between spatial coordinate identifiers and expression matrix column names")
  }
  if (length(common_cells) < ncol(data)) {
    warning(paste(
      "Spatial CellChat will use", length(common_cells), "of", ncol(data),
      "cells/spots with matched spatial coordinates"
    ))
  }
  spatial.df <- spatial.df[common_cells, , drop = FALSE]
  data <- data[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]

  use_spatial_factors <- !is.null(spatial.factors) && "spatial.factors" %in% names(formals(CellChat::computeCommunProb))
  if (!use_spatial_factors) {
    spatial.df <- offset_coordinates_by_library(spatial.df, interaction_length = opt$interaction_length)
  }
  spatial.mat <- as.matrix(spatial.df[, c("x", "y"), drop = FALSE])
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
future_max_size <- if (!is.null(opt$future_max_size_gb) && !is.na(opt$future_max_size_gb) && opt$future_max_size_gb > 0) {
  opt$future_max_size_gb * 1024^3
} else {
  Inf
}
options(future.globals.maxSize = future_max_size)
if (opt$nworkers <= 1) {
  future::plan("sequential")
} else {
  future::plan("multicore", workers = opt$nworkers)
}
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
if (opt$is_single_cell) {
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = opt$trim)
} else {
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

cell_group_sizes <- as.numeric(table(cellchat@idents))
names(cell_group_sizes) <- names(table(cellchat@idents))
network_height <- plot_dim(length(cell_group_sizes), base = 6.5, per = 0.18, min_size = 6.5, max_size = 10)
network_width <- network_height * 2
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_network.png")), width = network_width, height = network_height, units = "in", res = 300, bg = "white")
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = cell_group_sizes, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = cell_group_sizes, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
dev.off()
pdf(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_network.pdf")), width = network_width, height = network_height)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = cell_group_sizes, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = cell_group_sizes, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
dev.off()

png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_infoflow_bar.png")), width = 8.5, height = 6.5, units = "in", res = 300, bg = "white")
rankNet(cellchat, mode = "single", stacked = TRUE, do.stat = FALSE)
dev.off()

heatmap_size <- plot_dim(length(cell_group_sizes), base = 5.8, per = 0.28, min_size = 6, max_size = 12)
png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_heatmap.png")), width = heatmap_size * 2, height = heatmap_size, units = "in", res = 300, bg = "white")
h_count <- netVisual_heatmap(cellchat, measure = "count", title.name = "Number of interactions", font.size = 8, font.size.title = 10)
h_weight <- netVisual_heatmap(cellchat, measure = "weight", title.name = "Interaction strength", font.size = 8, font.size.title = 10)
print(h_count + h_weight)
dev.off()

cell_levels <- levels(cellchat@idents)
significant_lr <- net_lr
if ("prob" %in% names(significant_lr)) {
  significant_lr <- significant_lr[!is.na(significant_lr$prob) & significant_lr$prob > 0, , drop = FALSE]
}
if ("pval" %in% names(significant_lr)) {
  significant_lr <- significant_lr[!is.na(significant_lr$pval) & significant_lr$pval <= 0.05, , drop = FALSE]
}
if (all(c("source", "target") %in% names(significant_lr))) {
  significant_lr$cell_pair <- paste(significant_lr$source, significant_lr$target, sep = "|")
}

requested_cell_pairs <- parse_cell_pairs(opt$cell_pairs)
requested_focus_cells <- split_vec(opt$focus_cells)
requested_source_cells <- split_vec(opt$source_cells)
requested_target_cells <- split_vec(opt$target_cells)
selected_pair_keys <- NULL
selected_sources <- NULL
selected_targets <- NULL
selection_source <- "automatic"
scoped_lr <- significant_lr

if (!is.null(requested_cell_pairs)) {
  selection_source <- "cell_pairs"
  available_pair_keys <- unique(significant_lr$cell_pair)
  selected_pair_keys <- requested_cell_pairs[requested_cell_pairs %in% available_pair_keys]
  missing_pairs <- setdiff(requested_cell_pairs, selected_pair_keys)
  if (length(missing_pairs) > 0) {
    warning(paste("Requested CellChat cell pairs have no significant communication and will be skipped:", paste(missing_pairs, collapse = ", ")))
  }
  scoped_lr <- significant_lr[significant_lr$cell_pair %in% selected_pair_keys, , drop = FALSE]
} else if (!is.null(requested_focus_cells)) {
  selection_source <- "focus_cells"
  selected_focus_cells <- valid_named_values(opt$focus_cells, cell_levels, "focus cell")
  if (is.null(selected_focus_cells)) {
    scoped_lr <- significant_lr[0, , drop = FALSE]
  } else {
    selected_sources <- selected_focus_cells
    selected_targets <- selected_focus_cells
    scoped_lr <- significant_lr[
      significant_lr$source %in% selected_sources & significant_lr$target %in% selected_targets,
      , drop = FALSE
    ]
  }
} else if (!is.null(requested_source_cells) || !is.null(requested_target_cells)) {
  selection_source <- "source_target"
  selected_sources <- valid_named_values(opt$source_cells, cell_levels, "source cell")
  selected_targets <- valid_named_values(opt$target_cells, cell_levels, "target cell")
  invalid_source_scope <- !is.null(requested_source_cells) && is.null(selected_sources)
  invalid_target_scope <- !is.null(requested_target_cells) && is.null(selected_targets)
  if (invalid_source_scope || invalid_target_scope) {
    scoped_lr <- significant_lr[0, , drop = FALSE]
  } else {
    if (!is.null(selected_sources)) {
      scoped_lr <- scoped_lr[scoped_lr$source %in% selected_sources, , drop = FALSE]
    }
    if (!is.null(selected_targets)) {
      scoped_lr <- scoped_lr[scoped_lr$target %in% selected_targets, , drop = FALSE]
    }
  }
} else if (nrow(significant_lr) > 0) {
  pair_ranking <- significant_lr %>%
    dplyr::group_by(source, target, cell_pair) %>%
    dplyr::summarise(
      prob_sum = sum(prob, na.rm = TRUE),
      interaction_n = dplyr::n_distinct(interaction_name),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(prob_sum), dplyr::desc(interaction_n), source, target)
  selected_pair_keys <- head(pair_ranking$cell_pair, max(1, opt$top_cell_pairs))
  scoped_lr <- significant_lr[significant_lr$cell_pair %in% selected_pair_keys, , drop = FALSE]
}

if (nrow(scoped_lr) == 0) {
  warning("The selected CellChat cell scope contains no significant ligand-receptor communication; focused plots will be skipped")
}

scope_pathway_summary <- data.frame()
if (nrow(scoped_lr) > 0 && "pathway_name" %in% names(scoped_lr)) {
  scope_pathway_summary <- scoped_lr %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(
      interaction_n = dplyr::n(),
      lr_n = dplyr::n_distinct(interaction_name),
      source_n = dplyr::n_distinct(source),
      target_n = dplyr::n_distinct(target),
      prob_sum = sum(prob, na.rm = TRUE),
      prob_mean = mean(prob, na.rm = TRUE),
      prob_max = max(prob, na.rm = TRUE),
      pval_min = if ("pval" %in% names(scoped_lr)) min(pval, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(prob_sum), dplyr::desc(interaction_n), pathway_name)
}

requested_pathways <- split_vec(opt$pathways)
available_scope_pathways <- if (nrow(scope_pathway_summary) > 0) scope_pathway_summary$pathway_name else character(0)
selected_pathways <- select_pathways(opt$pathways, available_scope_pathways, scope_pathway_summary, opt$top_pathways)
pathway_selection_source <- if (!is.null(requested_pathways) && any(requested_pathways %in% selected_pathways)) "requested" else "automatic"

selected_lr_base <- scoped_lr
if (length(selected_pathways) > 0 && "pathway_name" %in% names(selected_lr_base)) {
  selected_lr_base <- selected_lr_base[selected_lr_base$pathway_name %in% selected_pathways, , drop = FALSE]
}

requested_lr_pairs <- parse_lr_pairs(opt$lr_pairs)
lr_selection_source <- "automatic"
selected_lr <- selected_lr_base
if (!is.null(requested_lr_pairs)) {
  requested_selected_lr <- filter_lr_pairs(selected_lr_base, requested_lr_pairs)
  if (nrow(requested_selected_lr) > 0) {
    selected_lr <- requested_selected_lr
    lr_selection_source <- "requested"
  } else {
    warning("Requested CellChat ligand-receptor pairs were not found in the selected cell scope; automatic top LR pairs will be used within the same scope")
    lr_selection_source <- "automatic_fallback"
  }
}

if (nrow(selected_lr) > 0) {
  if ("interaction_name" %in% names(selected_lr)) {
    selected_lr$interaction_key <- as.character(selected_lr$interaction_name)
  } else if ("interaction_name_2" %in% names(selected_lr)) {
    selected_lr$interaction_key <- as.character(selected_lr$interaction_name_2)
  } else {
    selected_lr$interaction_key <- as.character(selected_lr$lr_pair)
  }
  selected_lr <- selected_lr[!is.na(selected_lr$interaction_key) & nzchar(selected_lr$interaction_key), , drop = FALSE]
}

if (nrow(selected_lr) > 0) {
  lr_ranking <- selected_lr %>%
    dplyr::group_by(interaction_key) %>%
    dplyr::summarise(
      prob_sum = sum(prob, na.rm = TRUE),
      prob_max = max(prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(prob_sum), dplyr::desc(prob_max), interaction_key)
  coverage_keys <- character(0)
  if ("pathway_name" %in% names(selected_lr)) {
    for (pathway in selected_pathways) {
      pathway_lr <- selected_lr[selected_lr$pathway_name == pathway, , drop = FALSE]
      if (nrow(pathway_lr) > 0) {
        pathway_lr_rank <- pathway_lr %>%
          dplyr::group_by(interaction_key) %>%
          dplyr::summarise(prob_sum = sum(prob, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(dplyr::desc(prob_sum), interaction_key)
        coverage_keys <- c(coverage_keys, pathway_lr_rank$interaction_key[1])
      }
    }
  }
  if (!is.null(selected_pair_keys) && "cell_pair" %in% names(selected_lr)) {
    for (pair_key in selected_pair_keys) {
      pair_lr <- selected_lr[selected_lr$cell_pair == pair_key, , drop = FALSE]
      if (nrow(pair_lr) > 0) {
        pair_lr_rank <- pair_lr %>%
          dplyr::group_by(interaction_key) %>%
          dplyr::summarise(prob_sum = sum(prob, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(dplyr::desc(prob_sum), interaction_key)
        coverage_keys <- c(coverage_keys, pair_lr_rank$interaction_key[1])
      }
    }
  }
  coverage_keys <- unique(coverage_keys)
  ordered_interactions <- c(coverage_keys, setdiff(lr_ranking$interaction_key, coverage_keys))
  if (!is.null(opt$bubble_top_lr) && !is.na(opt$bubble_top_lr) && opt$bubble_top_lr > 0) {
    if (length(coverage_keys) > opt$bubble_top_lr) {
      warning("The requested CellChat scope needs more representative LR interactions than cellchat_bubble_top_lr allows; the highest-ranked representatives will be retained")
    }
    ordered_interactions <- head(ordered_interactions, opt$bubble_top_lr)
  }
  lr_ranking <- lr_ranking[match(ordered_interactions, lr_ranking$interaction_key), , drop = FALSE]
  selected_lr <- selected_lr[selected_lr$interaction_key %in% lr_ranking$interaction_key, , drop = FALSE]
  selected_lr$lr_rank <- match(selected_lr$interaction_key, lr_ranking$interaction_key)
  selected_lr <- selected_lr[order(selected_lr$lr_rank, -selected_lr$prob), , drop = FALSE]
} else {
  selected_lr$lr_rank <- integer(0)
}

selected_lr$selection_source <- rep(paste(selection_source, pathway_selection_source, lr_selection_source, sep = ":"), nrow(selected_lr))
selected_lr$pathway_rank <- if ("pathway_name" %in% names(selected_lr)) match(selected_lr$pathway_name, selected_pathways) else NA_integer_
selected_lr$plot_file <- rep(NA_character_, nrow(selected_lr))

selected_pathway_summary <- scope_pathway_summary[scope_pathway_summary$pathway_name %in% selected_pathways, , drop = FALSE]
if (nrow(selected_pathway_summary) > 0) {
  selected_pathway_summary$pathway_rank <- match(selected_pathway_summary$pathway_name, selected_pathways)
  selected_pathway_summary$selection_source <- paste(selection_source, pathway_selection_source, sep = ":")
  selected_pathway_summary <- selected_pathway_summary[order(selected_pathway_summary$pathway_rank), , drop = FALSE]
}

message(paste("CellChat focused selection source:", selection_source))
message(paste("Selected CellChat sources:", ifelse(is.null(selected_sources), "pair-specific or all", paste(selected_sources, collapse = ", "))))
message(paste("Selected CellChat targets:", ifelse(is.null(selected_targets), "pair-specific or all", paste(selected_targets, collapse = ", "))))
message(paste("Selected CellChat cell pairs:", ifelse(is.null(selected_pair_keys), "source-target cell scope", paste(selected_pair_keys, collapse = ", "))))
message(paste("Selected CellChat pathways:", ifelse(length(selected_pathways) == 0, "none", paste(selected_pathways, collapse = ", "))))
message(paste("Selected CellChat LR interactions:", ifelse(nrow(selected_lr) == 0, 0, dplyr::n_distinct(selected_lr$interaction_key))))

if (nrow(selected_lr) > 0 && !is.null(selected_pair_keys)) {
  plotted_pair_index <- 0
  for (pair_key in selected_pair_keys) {
    pair_lr <- selected_lr[selected_lr$cell_pair == pair_key, , drop = FALSE]
    if (nrow(pair_lr) == 0) {
      warning(paste("Skipping CellChat bubble plot because no selected LR remains for cell pair:", pair_key))
      next
    }
    plotted_pair_index <- plotted_pair_index + 1
    pair_parts <- unlist(strsplit(pair_key, "\\|"))
    bubble_basename <- if (plotted_pair_index == 1) {
      paste0(opt$sample_id, "_cellchat_bubble.png")
    } else {
      paste0(opt$sample_id, "_cellchat_bubble_", safe_label(pair_parts[1]), "_to_", safe_label(pair_parts[2]), ".png")
    }
    bubble_file <- file.path(opt$output_dir, bubble_basename)
    pair_lr_use <- make_pair_lr_use(pair_lr)
    bubble_height <- plot_dim(nrow(pair_lr_use), base = 3.5, per = 0.28, min_size = 4.5, max_size = 14)
    png(bubble_file, width = 8, height = bubble_height, units = "in", res = 300, bg = "white")
    p_bubble <- netVisual_bubble(
      cellchat,
      sources.use = pair_parts[1],
      targets.use = pair_parts[2],
      pairLR.use = pair_lr_use,
      remove.isolate = TRUE,
      sort.by.source = TRUE,
      sort.by.target = TRUE,
      angle.x = 45,
      font.size = 8,
      font.size.title = 10
    )
    print(p_bubble)
    dev.off()
    selected_lr$plot_file[selected_lr$cell_pair == pair_key] <- bubble_basename
  }
} else if (nrow(selected_lr) > 0) {
  bubble_basename <- paste0(opt$sample_id, "_cellchat_bubble.png")
  bubble_file <- file.path(opt$output_dir, bubble_basename)
  pair_lr_use <- make_pair_lr_use(selected_lr)
  bubble_width <- plot_dim(length(unique(c(selected_lr$source, selected_lr$target))), base = 6.5, per = 0.35, min_size = 7, max_size = 16)
  bubble_height <- plot_dim(nrow(pair_lr_use), base = 3.5, per = 0.28, min_size = 4.5, max_size = 14)
  png(bubble_file, width = bubble_width, height = bubble_height, units = "in", res = 300, bg = "white")
  p_bubble <- netVisual_bubble(
    cellchat,
    sources.use = selected_sources,
    targets.use = selected_targets,
    pairLR.use = pair_lr_use,
    remove.isolate = TRUE,
    sort.by.source = TRUE,
    sort.by.target = TRUE,
    angle.x = 45,
    font.size = 8,
    font.size.title = 10
  )
  print(p_bubble)
  dev.off()
  selected_lr$plot_file <- bubble_basename
}

selected_pair_plot <- NULL
if (!is.null(net_lr) && nrow(net_lr) > 0) {
  pair_candidates <- selected_lr
  if (!is.null(opt$pair_lr_use) && nzchar(opt$pair_lr_use)) {
    requested_pair <- normalise_lr_pair(opt$pair_lr_use)
    requested_pair_candidates <- filter_lr_pairs(pair_candidates, requested_pair)
    if (nrow(requested_pair_candidates) > 0) {
      pair_candidates <- requested_pair_candidates
    } else {
      warning(paste("Requested cellchat_pair_lr_use not found and automatic top LR pair will be used:", opt$pair_lr_use))
    }
  }
  if ("prob" %in% names(pair_candidates)) {
    pair_candidates <- pair_candidates[order(pair_candidates$prob, decreasing = TRUE), , drop = FALSE]
  }
  if (nrow(pair_candidates) > 0) {
    if ("interaction_name" %in% names(pair_candidates)) {
      selected_pair_plot <- as.character(pair_candidates$interaction_name[1])
    } else if ("interaction_name_2" %in% names(pair_candidates)) {
      selected_pair_plot <- as.character(pair_candidates$interaction_name_2[1])
    } else if ("lr_pair" %in% names(pair_candidates)) {
      selected_pair_plot <- as.character(pair_candidates$lr_pair[1])
    }
  }
}

if (!isTRUE(opt$is_single_cell) && !is.null(selected_pair_plot) && nzchar(selected_pair_plot)) {
  png(file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_spatial_lr.png")), width = 8, height = 6.5, units = "in", res = 300, bg = "white")
  p_spatial <- spatialFeaturePlot(cellchat, pairLR.use = selected_pair_plot, point.size = 0.8, do.binary = TRUE, cutoff = 0.05, enriched.only = FALSE, direction = 1, legend.size = 3, legend.text.size = 8)
  print(p_spatial)
  dev.off()
  message(paste("Spatial LR plot generated for:", selected_pair_plot))
}

requested_pathway_matches <- if (is.null(requested_pathways)) character(0) else intersect(requested_pathways, selected_pathways)
advanced_pathways <- if (length(requested_pathway_matches) > 0) selected_pathways else head(selected_pathways, 1)
if (isTRUE(opt$plot_advanced) && length(advanced_pathways) > 0) {
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  advanced_dir <- file.path(opt$output_dir, "advanced")
  dir.create(advanced_dir, showWarnings = FALSE, recursive = TRUE)
  for (pathway in advanced_pathways) {
    pathway_label <- safe_label(pathway)
    advanced_pair_keys <- if (!is.null(selected_pair_keys)) {
      if (selection_source == "automatic") head(selected_pair_keys, 1) else selected_pair_keys
    } else {
      NULL
    }
    if (!is.null(advanced_pair_keys)) {
      for (pair_key in advanced_pair_keys) {
        if (!any(selected_lr$cell_pair == pair_key & selected_lr$pathway_name == pathway)) {
          next
        }
        pair_parts <- unlist(strsplit(pair_key, "\\|"))
        pair_label <- paste0(safe_label(pair_parts[1]), "_to_", safe_label(pair_parts[2]))
        png(file.path(advanced_dir, paste0(pathway_label, "_", pair_label, "_aggregate_network.png")), width = 7, height = 6.5, units = "in", res = 300, bg = "white")
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", sources.use = pair_parts[1], targets.use = pair_parts[2], vertex.weight = cell_group_sizes, weight.scale = TRUE, vertex.label.cex = 0.8)
        dev.off()

        png(file.path(advanced_dir, paste0(pathway_label, "_", pair_label, "_lr_contribution.png")), width = 7.5, height = 5.5, units = "in", res = 300, bg = "white")
        p_contribution <- netAnalysis_contribution(cellchat, signaling = pathway, sources.use = pair_parts[1], targets.use = pair_parts[2], font.size = 8, font.size.title = 10)
        print(p_contribution)
        dev.off()
      }
    } else {
      png(file.path(advanced_dir, paste0(pathway_label, "_aggregate_network.png")), width = 7, height = 6.5, units = "in", res = 300, bg = "white")
      netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", sources.use = selected_sources, targets.use = selected_targets, vertex.weight = cell_group_sizes, weight.scale = TRUE, vertex.label.cex = 0.8)
      dev.off()

      png(file.path(advanced_dir, paste0(pathway_label, "_lr_contribution.png")), width = 7.5, height = 5.5, units = "in", res = 300, bg = "white")
      p_contribution <- netAnalysis_contribution(cellchat, signaling = pathway, sources.use = selected_sources, targets.use = selected_targets, font.size = 8, font.size.title = 10)
      print(p_contribution)
      dev.off()
    }

    png(file.path(advanced_dir, paste0(pathway_label, "_gene_expression.png")), width = 7.5, height = 5.5, units = "in", res = 300, bg = "white")
    p_gene <- plotGeneExpression(cellchat, signaling = pathway, type = "dot", color.use = c("#4C78A8", "#E07B73"))
    print(p_gene)
    dev.off()

    role_height <- plot_dim(length(cell_group_sizes), base = 4.5, per = 0.2, min_size = 5, max_size = 10)
    png(file.path(advanced_dir, paste0(pathway_label, "_signaling_role.png")), width = 8.5, height = role_height, units = "in", res = 300, bg = "white")
    netAnalysis_signalingRole_network(cellchat, signaling = pathway, slot.name = "netP", width = 8, height = role_height - 0.5, font.size = 8, font.size.title = 10)
    dev.off()
  }
}

saveRDS(cellchat, file.path(opt$output_dir, "cellchat.rds"))

write.csv(stats, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_stats.csv")), row.names = FALSE)
write.csv(net_lr, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_lr.csv")), row.names = FALSE)
write.csv(selected_lr, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_selected_lr.csv")), row.names = FALSE)
write.csv(selected_pathway_summary, file.path(opt$output_dir, paste0(opt$sample_id, "_cellchat_selected_pathway_summary.csv")), row.names = FALSE)
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
