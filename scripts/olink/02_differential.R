#!/usr/bin/env Rscript

if (!file.exists("README.md") || !file.exists("src/R/olink_analysis.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/olink_analysis.R")

cfg <- read_olink_config()
input_path <- "results/proteomics/olink/olink_preprocessed.rds"
if (!file.exists(input_path)) {
  stop("Run scripts/olink/01_qc.R first.")
}

long <- readRDS(input_path)
fit <- fit_olink_contrasts(long, cfg)

result_dir <- "results/proteomics/olink"
figure_dir <- "figures/exploratory/olink"
readr::write_csv(fit$results, file.path(result_dir, "differential_assays.csv"))

writeLines(
  c(
    paste0("formula: ~ ", cfg$model$fixed_effects),
    paste0("duplicateCorrelation consensus: ", fit$correlation)
  ),
  file.path(result_dir, "model_specification.txt")
)

p_volcano <- ggplot2::ggplot(
  fit$results,
  ggplot2::aes(x = .data$logFC, y = -log10(.data$P.Value))
) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::facet_wrap(~contrast, scales = "free_y") +
  ggplot2::labs(x = "log2 fold-change", y = "-log10(p-value)")

ggplot2::ggsave(
  file.path(figure_dir, "volcano.pdf"),
  plot = p_volcano,
  width = 10,
  height = 7
)

message("Olink differential analysis complete.")
