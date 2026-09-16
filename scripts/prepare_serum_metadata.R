#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Usage: Rscript scripts/prepare_serum_metadata.R [Metadata.csv]", call. = FALSE)
}

if (!file.exists("README.md") || !file.exists("src/R/serum_metadata.R")) {
  stop("Run this script from the repository root.")
}

source("src/R/serum_metadata.R")

metadata_path <- if (length(args) == 1L) args[[1L]] else "data/proteomics/Metadata.csv"
if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path)
}

out <- write_serum_metadata(metadata_path)

message("Mass-spec metadata: ", nrow(out$massspec), " assay-marked samples")
message("Olink metadata: ", nrow(out$olink), " assay-marked samples")
message("Wrote metadata/massspec_samples.csv and metadata/olink_samples.csv")
