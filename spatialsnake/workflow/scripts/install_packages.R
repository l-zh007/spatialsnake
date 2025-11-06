if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_packages <- c(
  "tidyverse", 
  "pheatmap", 
  "RColorBrewer", 
  "dplyr", 
  "optparse", 
  "ggplot2", 
  "ggsankey", 
  "cowplot")
bioc_packages <- c(
  "AnnotationDbi",
  "clusterProfiler", 
  "org.Hs.eg.db",
  "edgeR")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)}
}

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}