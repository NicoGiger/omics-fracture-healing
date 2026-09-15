#!/usr/bin/env Rscript

if (!file.exists("README.md") || !file.exists("src/R/proteomics_prolfqua.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/proteomics_prolfqua.R")

cfg <- read_proteomics_config()
lfq_path <- "results/proteomics/proteomics_lfq_preprocessed.rds"

if (!file.exists(lfq_path)) {
  stop("Preprocessed prolfqua object not found. Run scripts/01_proteomics_qc.R first.")
}

lfq <- load_prolfqua_bundle(lfq_path)
fit <- fit_proteomics_contrasts(lfq, cfg)

readr::write_csv(
  fit$results,
  "results/proteomics/differential_proteins.csv"
)

writeLines(
  c(
    paste0("mode: ", fit$mode),
    paste0("formula: ", fit$formula)
  ),
  "results/proteomics/model_specification.txt"
)

p_volcano <- fit$contrasts$get_Plotter()$volcano()$FDR
ggplot2::ggsave(
  "figures/exploratory/proteomics_volcano.pdf",
  plot = p_volcano,
  width = 9,
  height = 6
)

message("Differential analysis complete using ", fit$mode, ".")
message("Model: ", fit$formula)
