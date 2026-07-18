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

assert_version <- function(pkg, version) {
  if (!is_installed(pkg) || packageVersion(pkg) != version) {
    stop("Package '", pkg, "' expected version ", version, " but found ",
         if (is_installed(pkg)) as.character(packageVersion(pkg)) else "missing")
  }
}

github_matches <- function(pkg, version, sha) {
  if (!is_installed(pkg) || packageVersion(pkg) != version) return(FALSE)
  desc <- packageDescription(pkg)
  remote_sha <- desc$RemoteSha
  !is.null(remote_sha) && !is.na(remote_sha) && remote_sha == sha
}

install_github_pinned <- function(pkg, repo, version, sha, dependencies = r_dependencies) {
  if (github_matches(pkg, version, sha)) {
    message(pkg, " ", version, " at ", substr(sha, 1, 7), " is already installed")
    return(invisible(TRUE))
  }
  message("Installing ", pkg, " ", version, " from GitHub commit ", sha)
  remotes::install_github(repo, ref = sha, upgrade = "never", dependencies = dependencies)
  assert_version(pkg, version)
  assert_loadable(pkg)
}

install_cran_version <- function(pkg, version) {
  if (is_installed(pkg) && packageVersion(pkg) == version) {
    message(pkg, " ", version, " is already installed")
    return(invisible(TRUE))
  }
  message("Installing ", pkg, " ", version)
  remotes::install_version(pkg, version = version, upgrade = "never", dependencies = r_dependencies)
  assert_version(pkg, version)
  assert_loadable(pkg)
}

if (!is_installed("remotes")) install.packages("remotes", dependencies = NA)
if (!is_installed("BiocManager")) install.packages("BiocManager", dependencies = NA)

try({
  BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
}, silent = TRUE)
options(repos = BiocManager::repositories())

if (!is_installed("rhdf5")) {
  BiocManager::install("rhdf5", ask = FALSE, update = FALSE)
}
assert_loadable("rhdf5")

install_github_pinned(
  pkg = "ggsankey",
  repo = "davidsjoberg/ggsankey",
  version = "0.0.99999",
  sha = "b675d0d5144b1b5758d3b2b41e86ceee66a1e071"
)

install_github_pinned(
  pkg = "schard",
  repo = "cellgeni/schard",
  version = "1.1.0",
  sha = "f9c30c25a192175a1e6523e32f0e2047de1e8685",
  dependencies = FALSE
)

install_github_pinned(
  pkg = "spacexr",
  repo = "dmcable/spacexr",
  version = "2.2.1",
  sha = "9f5dc33c8060f946c6072a138b70e189636e1435"
)

install_cran_version("NMF", "0.26")
install_cran_version("circlize", "0.4.15")

if (!is_installed("presto")) {
  remotes::install_github("immunogenomics/presto", upgrade = "never", dependencies = r_dependencies)
}
assert_loadable("presto")

install_github_pinned(
  pkg = "CellChat",
  repo = "jinworks/CellChat",
  version = "2.2.0.9001",
  sha = "75253cd0c9e68410e6e721a6d3a0419a1d7e358f"
)

required <- c("ggsankey", "schard", "spacexr", "CellChat", "NMF", "circlize")
miss <- required[!vapply(required, is_installed, logical(1))]
if (length(miss) > 0) {
  stop("Missing packages after install: ", paste(miss, collapse = ", "))
}

for (pkg in required) {
  assert_loadable(pkg)
}

message("All conda GitHub and version-pinned R packages are installed and loadable.")
