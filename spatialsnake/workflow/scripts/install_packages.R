#!/usr/bin/env Rscript

Sys.setenv(C_INCLUDE_PATH = "")
Sys.setenv(CPLUS_INCLUDE_PATH = "")
Sys.setenv(CPATH = "")

r <- getOption("repos")
r["CRAN"] <- "https://cloud.r-project.org"
options(repos = r)
options(timeout = 600000)

.libPaths(unique(c(.Library, .libPaths())))

r_dependencies <- c("Depends", "Imports", "LinkingTo")

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

if (!is_installed("BiocManager")) install.packages("BiocManager", dependencies = NA)
try({
  BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
}, silent = TRUE)
options(repos = BiocManager::repositories())
BiocManager::install(c("rhdf5"), ask = FALSE, update = FALSE)

required_schard <- "1.1.0"
schard_ref <- "f9c30c25a192175a1e6523e32f0e2047de1e8685"
if (!is_installed("schard") || packageVersion("schard") != required_schard) {
  message("Installing tested schard ", required_schard, " from GitHub")
  remotes::install_github(
    paste0("cellgeni/schard@", schard_ref),
    upgrade = "never",
    dependencies = FALSE,
    force = TRUE
  )
}
if (packageVersion("schard") != required_schard) {
  stop("schard ", required_schard, " is required after installation.")
}
spacexr_ref <- "9f5dc33c8060f946c6072a138b70e189636e1435"
if (!is_installed("spacexr") || packageVersion("spacexr") != "2.2.1" ||
    !identical(packageDescription("spacexr")$RemoteSha, spacexr_ref)) {
  message("Installing tested spacexr 2.2.1 from GitHub")
  remotes::install_github(paste0("dmcable/spacexr@", spacexr_ref), upgrade = "never", dependencies = r_dependencies, build_vignettes = FALSE, force = TRUE)
}

required_nmf <- "0.26"
required_circlize <- "0.4.15"

if (is_installed("NMF")) {
  if (packageVersion("NMF") != required_nmf) {
    message("Downgrading NMF to ", required_nmf)
    remotes::install_version("NMF", version = required_nmf, upgrade = "never", dependencies = r_dependencies)
  }
} else {
  remotes::install_version("NMF", version = required_nmf, upgrade = "never", dependencies = r_dependencies)
}

if (is_installed("circlize")) {
  if (packageVersion("circlize") != required_circlize) {
    message("Downgrading circlize to ", required_circlize)
    remotes::install_version("circlize", version = required_circlize, upgrade = "never", dependencies = r_dependencies)
  }
} else {
  remotes::install_version("circlize", version = required_circlize, upgrade = "never", dependencies = r_dependencies)
}

if (!is_installed("ggpubr")) {
  stop(
    "CellChat requires ggpubr. Install the tested Conda package first with: ",
    "conda install --override-channels -c conda-forge -c bioconda r-ggpubr -y"
  )
}
assert_loadable("ggpubr")

cellchat_ref <- "75253cd0c9e68410e6e721a6d3a0419a1d7e358f"
if (!is_installed("CellChat") || packageVersion("CellChat") != "2.2.0.9001" ||
    !identical(packageDescription("CellChat")$RemoteSha, cellchat_ref)) {
  message("Installing tested CellChat 2.2.0.9001 from GitHub")
  if (!is_installed("presto")) {
    remotes::install_github("immunogenomics/presto", upgrade = "never", dependencies = r_dependencies)
  }
  remotes::install_github(paste0("jinworks/CellChat@", cellchat_ref), upgrade = "never", dependencies = r_dependencies, force = TRUE)
}

required <- c("schard", "spacexr", "CellChat", "ggpubr", "NMF", "circlize")
miss <- required[!vapply(required, is_installed, logical(1))]
if (length(miss) > 0) {
  stop("Missing packages after install: ", paste(miss, collapse = ", "))
}

for (pkg in required) {
  assert_loadable(pkg)
}

message("All GitHub and version‑pinned packages are installed and loadable.")
