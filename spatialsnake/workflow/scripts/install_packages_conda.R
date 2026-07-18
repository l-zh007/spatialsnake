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
seurat_dependencies <- c("Depends", "Imports", "LinkingTo")

normalize_version <- function(version) {
  gsub("-", ".", version, fixed = TRUE)
}

install_cran_version <- function(pkg, version) {
  if (is_installed(pkg) && as.character(packageVersion(pkg)) == normalize_version(version)) {
    message(pkg, " ", version, " is already installed")
    return(invisible(TRUE))
  }
  message("Installing ", pkg, " ", version)
  remotes::install_version(pkg, version = version, upgrade = "never", dependencies = seurat_dependencies)
}

if (is_installed("Seurat")) {
  if (packageVersion("Seurat") != required_seurat) {
    message("Installing Seurat ", required_seurat)
    remotes::install_version("Seurat", version = required_seurat, upgrade = "never", dependencies = seurat_dependencies)
  }
} else {
  remotes::install_version("Seurat", version = required_seurat, upgrade = "never", dependencies = seurat_dependencies)
}

install_cran_version("spam", "2.11-3")
install_cran_version("spatstat.random", "3.4-5")
install_cran_version("spatstat.explore", "3.8-0")
install_cran_version("leidenbase", "0.1.32")

assert_loadable("Seurat")
assert_loadable("leidenbase")

message("Conda-only R packages are installed and loadable.")
