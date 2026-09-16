#!/usr/bin/env Rscript

if (!file.exists("README.md") || !file.exists("src/R/massspec_prolfqua.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/massspec_prolfqua.R")

cfg <- read_ms_config()
lfq_path <- "results/proteomics/massspec/ms_lfq_preprocessed.rds"
if (!file.exists(lfq_path)) {
  stop("Run scripts/01_ms_qc.R first.")
}

lfq <- load_ms_lfq(lfq_path)
bio_lfq <- ms_biological_lfq(lfq, cfg)
fit <- fit_ms_contrasts(bio_lfq, cfg)

result_dir <- "results/proteomics/massspec"
figure_dir <- "figures/exploratory/massspec"

readr::write_csv(fit$results, file.path(result_dir, "differential_proteins.csv"))
writeLines(
  c(
    paste0("mode: ", fit$mode),
    paste0("formula: ", fit$formula)
  ),
  file.path(result_dir, "model_specification.txt")
)

p_volcano <- fit$contrasts$get_Plotter()$volcano()$FDR
ggplot2::ggsave(
  file.path(figure_dir, "volcano.pdf"),
  plot = p_volcano,
  width = 9,
  height = 6
)

message("MS differential analysis complete using ", fit$mode, ".")
message("Model: ", fit$formula)
