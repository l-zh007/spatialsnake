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

if (!is_installed("BiocManager")) install.packages("BiocManager", dependencies = NA)
try({
  BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
}, silent = TRUE)
options(repos = BiocManager::repositories())
BiocManager::install(c("rhdf5"), ask = FALSE, update = FALSE)

if (!is_installed("ggsankey")) {
  message("Installing ggsankey from GitHub")
  remotes::install_github("davidsjoberg/ggsankey", upgrade = "never", dependencies = TRUE)
}
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
if (!is_installed("spacexr")) {
  message("Installing spacexr from GitHub")
  remotes::install_github("dmcable/spacexr", upgrade = "never", dependencies = TRUE, build_vignettes = FALSE)
}

required_nmf <- "0.26"
required_circlize <- "0.4.15"

if (is_installed("NMF")) {
  if (packageVersion("NMF") != required_nmf) {
    message("Downgrading NMF to ", required_nmf)
    remotes::install_version("NMF", version = required_nmf, upgrade = "never", dependencies = TRUE)
  }
} else {
  remotes::install_version("NMF", version = required_nmf, upgrade = "never", dependencies = TRUE)
}

if (is_installed("circlize")) {
  if (packageVersion("circlize") != required_circlize) {
    message("Downgrading circlize to ", required_circlize)
    remotes::install_version("circlize", version = required_circlize, upgrade = "never", dependencies = TRUE)
  }
} else {
  remotes::install_version("circlize", version = required_circlize, upgrade = "never", dependencies = TRUE)
}

if (!is_installed("CellChat")) {
  message("Installing CellChat from GitHub")
  if (!is_installed("presto")) {
    remotes::install_github("immunogenomics/presto", upgrade = "never", dependencies = TRUE)
  }
  remotes::install_github("jinworks/CellChat", upgrade = "never", dependencies = TRUE)
}

required <- c("ggsankey", "schard", "spacexr", "CellChat", "NMF", "circlize")
miss <- required[!vapply(required, is_installed, logical(1))]
if (length(miss) > 0) {
  stop("Missing packages after install: ", paste(miss, collapse = ", "))
}

for (pkg in required) {
  assert_loadable(pkg)
}

message("All GitHub and version‑pinned packages are installed and loadable.")
