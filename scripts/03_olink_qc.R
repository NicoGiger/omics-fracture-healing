#!/usr/bin/env Rscript

if (!file.exists("README.md") || !file.exists("src/R/olink_analysis.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/olink_analysis.R")

cfg <- read_olink_config()
paths <- read_olink_paths()
parsed <- read_olink_quant(paths$quant, paths$metadata, cfg)
long <- prepare_olink_data(parsed$long, cfg)
qc <- olink_qc_tables(long)

result_dir <- "results/proteomics/olink"
figure_dir <- "figures/exploratory/olink"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(parsed$annotations, file.path(result_dir, "assay_annotations.csv"))
readr::write_csv(qc$sample, file.path(result_dir, "qc_sample_summary.csv"))
readr::write_csv(qc$assay, file.path(result_dir, "qc_assay_missingness.csv"))
readr::write_csv(qc$design, file.path(result_dir, "sample_design_counts.csv"))

writeLines(
  olink_design_columns(long, cfg),
  file.path(result_dir, "model_matrix_columns.txt")
)

saveRDS(long, file.path(result_dir, "olink_preprocessed.rds"))

p_missing <- ggplot2::ggplot(
  qc$assay,
  ggplot2::aes(x = stats::reorder(.data$assay, .data$missing_fraction), y = .data$missing_fraction)
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = "Assay", y = "Missing fraction")

ggplot2::ggsave(
  file.path(figure_dir, "assay_missingness.pdf"),
  plot = p_missing,
  width = 7,
  height = 8
)

p_density <- long |>
  dplyr::filter(!is.na(.data$log2_concentration)) |>
  ggplot2::ggplot(ggplot2::aes(x = .data$log2_concentration, group = .data$sample_id)) +
  ggplot2::geom_density(alpha = 0.15) +
  ggplot2::labs(x = "log2 concentration (pg/mL)", y = "Density")

ggplot2::ggsave(
  file.path(figure_dir, "concentration_density.pdf"),
  plot = p_density,
  width = 8,
  height = 5
)

message("Olink QC complete. Results: ", result_dir)
