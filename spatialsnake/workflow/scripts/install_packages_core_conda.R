#!/usr/bin/env Rscript

.libPaths(unique(c(.Library, .libPaths())))
required <- c("AnnotationDbi", "dplyr", "clusterProfiler", "org.Hs.eg.db",
              "org.Mm.eg.db", "optparse", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing core R packages: ", paste(missing, collapse = ", "))
message("Core R enrichment packages are installed and loadable.")
