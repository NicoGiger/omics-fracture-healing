#!/usr/bin/env Rscript

if (!file.exists("README.md") || !file.exists("src/R/massspec_prolfqua.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/massspec_prolfqua.R")

cfg <- read_ms_config()
paths <- read_ms_paths()
input <- read_ms_wide(paths$abundance, paths$metadata, cfg)

# Keep the analysis path close to prolfqua's standard workflow:
# AnalysisConfiguration -> setup_analysis -> LFQData -> Transformer -> Plotter.
lfq <- make_ms_lfq(input$data, cfg)
lfq <- preprocess_ms_lfq(lfq, cfg)
qc <- ms_qc_tables(lfq, cfg)
bio_lfq <- ms_biological_lfq(lfq, cfg)

result_dir <- "results/proteomics/massspec"
figure_dir <- "figures/exploratory/massspec"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(qc$sample, file.path(result_dir, "qc_sample_missingness.csv"))
readr::write_csv(qc$protein, file.path(result_dir, "qc_protein_missingness.csv"))
readr::write_csv(
  input$metadata,
  file.path(result_dir, "sample_metadata_used.csv")
)

writeLines(
  ms_design_columns(bio_lfq, cfg),
  file.path(result_dir, "model_matrix_columns.txt")
)

save_ms_lfq(lfq, file.path(result_dir, "ms_lfq_preprocessed.rds"))

p_density <- lfq$get_Plotter()$intensity_distribution_density()
ggplot2::ggsave(
  file.path(figure_dir, "intensity_density.pdf"),
  plot = p_density,
  width = 8,
  height = 5
)

message("MS response: ", lfq$response())
message("Repeated biological animals detected: ", has_repeated_ms_animals(bio_lfq, cfg))
message("MS QC complete. Results: ", result_dir)
