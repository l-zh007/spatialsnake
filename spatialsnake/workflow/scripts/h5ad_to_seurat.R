#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Default to 'st' if type is not provided, for robustness
data_type <- "st"
if (length(args) >= 3) {
  data_type <- args[3]
} else if (length(args) < 2) {
  stop("Usage: Rscript h5ad_to_seurat.R <input_h5ad> <output_rds> [type]")
}

input_h5ad <- args[1]
output_rds <- args[2]

if (!requireNamespace("schard", quietly = TRUE)) {
  stop("Package 'schard' is not installed. Please install it first.")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Package 'Seurat' is not installed. Please install it first.")
}

library(schard)
library(Seurat)

message(paste("Converting", input_h5ad, "to", output_rds))
message(paste("Data type:", data_type))

if (data_type == "sc") {
  message("Using schard::h5ad2seurat for Single Cell data...")
  obj <- schard::h5ad2seurat(input_h5ad, use.raw = TRUE)
  
} else {
  message("Using schard::h5ad2seurat_spatial for Spatial Transcriptomics data...")
  
  obj <- tryCatch({
    message("Attempting standard spatial conversion...")
    schard::h5ad2seurat_spatial(input_h5ad, use.raw = FALSE)
    
  }, error = function(e) {
    message("Standard spatial conversion failed or encountered error.")
    message(paste("Error message:", e$message))
    message("Attempting multi-sample conversion (simplify = FALSE) with merge loop...")
    
    visl <- schard::h5ad2seurat_spatial(input_h5ad, use.raw = FALSE, simplify = FALSE)
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
}

saveRDS(obj, file = output_rds)
message("Conversion complete.")
