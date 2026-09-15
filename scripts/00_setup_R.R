#!/usr/bin/env Rscript

if (!file.exists("README.md") || !dir.exists("src/R")) {
  stop("Run this script from the repository root.")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE, restart = FALSE)
}

renv::install(c(
  "yaml",
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2"
))

# Pin the stable release rather than the moving development branch.
renv::install("fgcz/prolfqua@v.1.5.0")
renv::snapshot()

message("R environment initialized and snapshotted.")
