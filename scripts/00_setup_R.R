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
  "readxl",
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2"
))

# prolfqua is used only for the DIA-NN / mass-spec branch.
renv::install("fgcz/prolfqua@v.1.5.0")
# Olink quantified-concentration modelling uses limma directly.
renv::install("bioc::limma")
renv::snapshot()

message("R environment initialized and snapshotted.")
