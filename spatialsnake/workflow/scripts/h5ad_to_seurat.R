#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    paste(
      "Usage: Rscript h5ad_to_seurat.R",
      "<input_h5ad> <output_rds> <st|sc> <counts|data>",
      "<spatial|generic> <minimum_schard_version>"
    )
  )
}

input_h5ad <- args[[1]]
output_rds <- args[[2]]
data_type <- tolower(args[[3]])
target_layer <- tolower(args[[4]])
conversion_mode <- tolower(args[[5]])
minimum_schard_version <- args[[6]]

if (!file.exists(input_h5ad)) {
  stop(sprintf("Input H5AD does not exist: %s", input_h5ad))
}
if (!data_type %in% c("st", "sc")) {
  stop("data_type must be 'st' or 'sc'.")
}
if (!target_layer %in% c("counts", "data")) {
  stop("target_layer must be 'counts' or 'data'.")
}
if (!conversion_mode %in% c("spatial", "generic")) {
  stop("conversion_mode must be 'spatial' or 'generic'.")
}

for (package in c("schard", "Seurat", "SeuratObject", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Required R package '%s' is not installed.", package))
  }
}

installed_schard <- utils::packageVersion("schard")
if (installed_schard < numeric_version(minimum_schard_version)) {
  stop(
    sprintf(
      paste0(
        "Schard >= %s is required, but version %s is installed. ",
        "Run 'spatialsnake install-packages' to install the tested release."
      ),
      minimum_schard_version,
      as.character(installed_schard)
    )
  )
}

message(sprintf("Schard version: %s", as.character(installed_schard)))
message(sprintf("Input data type: %s", data_type))
message(sprintf("Conversion mode: %s", conversion_mode))
message(sprintf("Target Seurat layer: %s", target_layer))

if (conversion_mode == "spatial") {
  obj <- schard::h5ad2seurat_spatial(
    input_h5ad,
    use.raw = FALSE,
    simplify = TRUE
  )
} else {
  obj <- schard::h5ad2seurat(
    input_h5ad,
    use.raw = FALSE,
    load.obsm = TRUE,
    load.X = TRUE
  )
}

if (!inherits(obj, "Seurat")) {
  stop("Schard did not return a Seurat object.")
}
if (ncol(obj) == 0L || nrow(obj) == 0L) {
  stop("Schard returned an empty Seurat object.")
}

assay_name <- SeuratObject::DefaultAssay(obj)
old_assay <- obj[[assay_name]]
if (inherits(old_assay, "Assay5")) {
  obj <- SeuratObject::JoinLayers(obj, assay = assay_name)
  old_assay <- obj[[assay_name]]
}
source_layers <- SeuratObject::Layers(old_assay)
source_layer <- if (any(source_layers == "counts")) "counts" else source_layers[[1]]
values <- SeuratObject::LayerData(obj, assay = assay_name, layer = source_layer)
if (!inherits(values, "sparseMatrix")) {
  values <- Matrix::Matrix(values, sparse = TRUE)
}
values <- methods::as(methods::as(values, "generalMatrix"), "CsparseMatrix")
feature_metadata <- old_assay[[]]
assay_key <- SeuratObject::Key(old_assay)

if (target_layer == "data") {
  data_assay <- SeuratObject::CreateAssay5Object(data = values)
  replacement_assay <- data_assay
} else {
  counts_assay <- SeuratObject::CreateAssay5Object(counts = values)
  replacement_assay <- counts_assay
}
replacement_assay <- SeuratObject::AddMetaData(
  replacement_assay,
  metadata = feature_metadata
)
SeuratObject::Key(replacement_assay) <- assay_key
obj[[assay_name]] <- replacement_assay

if (conversion_mode == "spatial" && length(SeuratObject::Images(obj)) == 0L) {
  stop("Spatial Schard conversion completed without creating a Seurat image.")
}

assay_layers <- SeuratObject::Layers(obj[[assay_name]], search = NA)
if (target_layer == "counts") {
  if (!any(startsWith(assay_layers, "counts"))) {
    stop("Count-matrix conversion completed without a counts layer.")
  }
  if (any(startsWith(assay_layers, "data"))) {
    stop("Count-matrix conversion unexpectedly created a duplicate data layer.")
  }
} else {
  if (!any(startsWith(assay_layers, "data"))) {
    stop("Normalized-matrix conversion completed without a data layer.")
  }
  if (any(startsWith(assay_layers, "counts"))) {
    stop("Normalized-matrix conversion unexpectedly created a counts layer.")
  }
}

saveRDS(obj, file = output_rds)
validated <- readRDS(output_rds)
if (!inherits(validated, "Seurat") || ncol(validated) != ncol(obj) || nrow(validated) != nrow(obj)) {
  stop("RDS validation failed after writing the converted Seurat object.")
}
message(sprintf("Validated Seurat RDS: %d observations, %d features", ncol(obj), nrow(obj)))
