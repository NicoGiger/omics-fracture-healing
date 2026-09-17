#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop(
    "Usage: Rscript scripts/spatial/01_import_qc.R [spatial_data_root]",
    call. = FALSE
  )
}

if (!file.exists("README.md") || !file.exists("src/R/spatial_utils.R")) {
  stop("Run this script from the repository root.", call. = FALSE)
}

required_packages <- c("Seurat", "SeuratObject", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the project renv before running this script.",
    call. = FALSE
  )
}

source("src/R/spatial_utils.R")

# Keep QC choices visible here rather than hiding them in another config file.
min_features <- 100L
max_percent_mt <- 10

data_root <- if (length(args) == 1L) {
  args[[1L]]
} else {
  Sys.getenv("SPATIAL_DATA_ROOT", unset = "data/spatial/raw")
}
if (!dir.exists(data_root)) {
  stop(
    "Spatial data root not found: ", data_root,
    ". Pass it as the first argument or set SPATIAL_DATA_ROOT.",
    call. = FALSE
  )
}
data_root <- normalizePath(data_root, mustWork = TRUE)

samples <- read_spatial_samples("metadata/spatial_samples.csv")

object_dir <- "data/spatial/processed"
result_dir <- "results/spatial/qc"
figure_dir <- "figures/exploratory/spatial/qc"

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sample_qc <- vector("list", nrow(samples))
annotation_coverage <- vector("list", nrow(samples))

for (i in seq_len(nrow(samples))) {
  capture <- samples[i, , drop = FALSE]
  message(
    "[", i, "/", nrow(samples), "] ",
    capture$capture_id, " (", capture$time_point, ")"
  )

  out <- process_spatial_capture(
    capture = capture,
    data_root = data_root,
    min_features = min_features,
    max_percent_mt = max_percent_mt
  )

  saveRDS(
    out$object,
    file.path(object_dir, paste0(capture$capture_id, ".rds"))
  )

  save_spatial_qc_plots(
    out$object,
    capture_id = capture$capture_id,
    output_dir = figure_dir
  )

  sample_qc[[i]] <- out$sample_qc
  annotation_coverage[[i]] <- out$annotation_coverage
}

sample_qc <- do.call(rbind, sample_qc)
annotation_coverage <- do.call(rbind, annotation_coverage)

utils::write.csv(
  sample_qc,
  file.path(result_dir, "sample_qc.csv"),
  row.names = FALSE
)
utils::write.csv(
  annotation_coverage,
  file.path(result_dir, "annotation_coverage.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "session_info.txt")
)

message("Spatial import and QC complete.")
message("Processed objects: ", object_dir)
message("QC tables: ", result_dir)
message("QC figures: ", figure_dir)
