#!/usr/bin/env Rscript

if (!file.exists("README.md") || !file.exists("src/R/proteomics_prolfqua.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/proteomics_prolfqua.R")

cfg <- read_proteomics_config()
paths <- read_proteomics_paths()

input <- read_proteomics_wide(
  abundance_path = paths$abundance,
  metadata_path = paths$metadata,
  cfg = cfg
)

lfq <- make_prolfqua_data(input$data, cfg)
lfq <- preprocess_prolfqua(lfq, cfg)
qc <- proteomics_qc_tables(lfq, cfg)

dir.create("results/proteomics", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/exploratory", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(qc$sample, "results/proteomics/qc_sample_missingness.csv")
readr::write_csv(qc$protein, "results/proteomics/qc_protein_missingness.csv")

writeLines(
  proteomics_design_columns(lfq, cfg),
  "results/proteomics/model_matrix_columns.txt"
)

save_prolfqua_bundle(
  lfq,
  "results/proteomics/proteomics_lfq_preprocessed.rds"
)

p_density <- lfq$get_Plotter()$intensity_distribution_density()
ggplot2::ggsave(
  "figures/exploratory/proteomics_intensity_density.pdf",
  plot = p_density,
  width = 8,
  height = 5
)

message("Response: ", lfq$response())
message("Repeated animals detected: ", has_repeated_animals(lfq, cfg))
message("QC outputs written to results/proteomics/ and figures/exploratory/.")
