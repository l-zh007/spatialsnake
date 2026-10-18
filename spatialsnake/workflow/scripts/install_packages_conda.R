#!/usr/bin/env Rscript

Sys.setenv(C_INCLUDE_PATH = "")
Sys.setenv(CPLUS_INCLUDE_PATH = "")
Sys.setenv(CPATH = "")

r <- getOption("repos")
r["CRAN"] <- "https://cloud.r-project.org"
options(repos = r)
options(timeout = 600000)

.libPaths(unique(c(.Library, .libPaths())))

is_installed <- function(pkg) {
  pkg %in% rownames(installed.packages(lib.loc = .libPaths()))
}

assert_loadable <- function(pkg) {
  tryCatch(
    {
      loadNamespace(pkg)
      invisible(TRUE)
    },
    error = function(e) {
      stop("Package '", pkg, "' is installed but cannot be loaded: ", conditionMessage(e))
    }
  )
}

if (!is_installed("remotes")) install.packages("remotes", dependencies = NA)

required_seurat <- "5.5.0"
conda_seurat <- "5.5.1"
seurat_dependencies <- c("Depends", "Imports", "LinkingTo")

if (!is_installed("Seurat")) {
  remotes::install_version("Seurat", version = required_seurat, upgrade = "never", dependencies = seurat_dependencies)
} else if (!as.character(packageVersion("Seurat")) %in% c(required_seurat, conda_seurat)) {
  stop("Seurat ", required_seurat, " or ", conda_seurat, " is required")
}

for (pkg in c("Seurat", "spam", "spatstat.random", "spatstat.explore", "leidenbase")) {
  assert_loadable(pkg)
}

message("Conda-only R packages are installed and loadable.")
