#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L || length(args) > 2L) {
  stop(
    paste0(
      "Usage: Rscript scripts/massspec/01_extract_diann_quant.R <report.tsv> [output_dir]\n",
      "Example: Rscript scripts/massspec/01_extract_diann_quant.R /path/to/report.tsv data/proteomics/massspec/diann_extract"
    ),
    call. = FALSE
  )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required. Run Rscript scripts/setup_R.R first.",
    call. = FALSE
  )
}

input_file <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- if (length(args) == 2L) args[[2L]] else "data/proteomics/massspec/diann_extract"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

header <- names(
  data.table::fread(
    input_file,
    nrows = 0L,
    check.names = FALSE,
    showProgress = FALSE
  )
)

writeLines(header, file.path(output_dir, "diann_columns.txt"))

candidate_columns <- c(
  "Run",
  "Protein.Group",
  "Genes",
  "Protein.Names",
  "PG.MaxLFQ",
  "PG.Normalised",
  "PG.Quantity",
  "PG.Q.Value",
  "Q.Value",
  "Global.Q.Value",
  "Proteotypic",
  "Precursor.Id"
)

required_columns <- c("Run", "Protein.Group")
missing_required <- setdiff(required_columns, header)
if (length(missing_required) > 0L) {
  stop(
    "DIA-NN report is missing required columns: ",
    paste(missing_required, collapse = ", "),
    call. = FALSE
  )
}

available_columns <- intersect(candidate_columns, header)
quant_columns <- intersect(
  c("PG.MaxLFQ", "PG.Normalised", "PG.Quantity"),
  available_columns
)

if (length(quant_columns) == 0L) {
  warning(
    "None of PG.MaxLFQ, PG.Normalised, or PG.Quantity was found. ",
    "The column list will still be written for inspection."
  )
}

message("Reading selected DIA-NN columns: ", paste(available_columns, collapse = ", "))
dt <- data.table::fread(
  input_file,
  select = available_columns,
  na.strings = c("", "NA", "NaN"),
  check.names = FALSE,
  showProgress = TRUE
)

first_non_missing <- function(x) {
  idx <- which(!is.na(x))
  if (length(idx) > 0L) {
    return(x[[idx[[1L]]]])
  }
  x[NA_integer_][[1L]]
}

carry_columns <- intersect(
  c(
    "Genes",
    "Protein.Names",
    "PG.MaxLFQ",
    "PG.Normalised",
    "PG.Quantity",
    "PG.Q.Value"
  ),
  names(dt)
)

protein_quant <- if (length(carry_columns) > 0L) {
  dt[
    , lapply(.SD, first_non_missing),
    by = .(Run, Protein.Group),
    .SDcols = carry_columns
  ]
} else {
  unique(dt[, .(Run, Protein.Group)])
}

row_counts <- dt[, .(n_report_rows = .N), by = .(Run, Protein.Group)]
protein_quant <- merge(
  protein_quant,
  row_counts,
  by = c("Run", "Protein.Group"),
  all.x = TRUE,
  sort = FALSE
)

if ("Precursor.Id" %in% names(dt)) {
  precursor_counts <- dt[
    , .(n_precursors = data.table::uniqueN(Precursor.Id, na.rm = TRUE)),
    by = .(Run, Protein.Group)
  ]
  protein_quant <- merge(
    protein_quant,
    precursor_counts,
    by = c("Run", "Protein.Group"),
    all.x = TRUE,
    sort = FALSE
  )
}

if ("Q.Value" %in% names(dt)) {
  q_counts <- dt[
    , .(n_q01 = sum(Q.Value < 0.01, na.rm = TRUE)),
    by = .(Run, Protein.Group)
  ]
  protein_quant <- merge(
    protein_quant,
    q_counts,
    by = c("Run", "Protein.Group"),
    all.x = TRUE,
    sort = FALSE
  )
}

if ("Global.Q.Value" %in% names(dt)) {
  global_q_counts <- dt[
    , .(n_global_q01 = sum(Global.Q.Value < 0.01, na.rm = TRUE)),
    by = .(Run, Protein.Group)
  ]
  protein_quant <- merge(
    protein_quant,
    global_q_counts,
    by = c("Run", "Protein.Group"),
    all.x = TRUE,
    sort = FALSE
  )
}

if ("Proteotypic" %in% names(dt)) {
  proteotypic_counts <- dt[
    , .(n_proteotypic = sum(Proteotypic == 1, na.rm = TRUE)),
    by = .(Run, Protein.Group)
  ]
  protein_quant <- merge(
    protein_quant,
    proteotypic_counts,
    by = c("Run", "Protein.Group"),
    all.x = TRUE,
    sort = FALSE
  )
}

if (all(c("Q.Value", "Proteotypic") %in% names(dt))) {
  proteotypic_q_counts <- dt[
    , .(
      n_proteotypic_q01 = sum(
        Proteotypic == 1 & Q.Value < 0.01,
        na.rm = TRUE
      )
    ),
    by = .(Run, Protein.Group)
  ]
  protein_quant <- merge(
    protein_quant,
    proteotypic_q_counts,
    by = c("Run", "Protein.Group"),
    all.x = TRUE,
    sort = FALSE
  )
}

# Protein-level DIA-NN fields should normally be constant across the precursor
# rows belonging to one Run x Protein.Group. Record violations instead of
# silently assuming that this is always true.
consistency_columns <- intersect(
  c(
    "Genes",
    "Protein.Names",
    "PG.MaxLFQ",
    "PG.Normalised",
    "PG.Quantity",
    "PG.Q.Value"
  ),
  names(dt)
)

consistency <- data.table::rbindlist(
  lapply(consistency_columns, function(column_name) {
    tmp <- dt[
      , .(
        n_unique = data.table::uniqueN(
          get(column_name),
          na.rm = TRUE
        )
      ),
      by = .(Run, Protein.Group)
    ]

    data.table::data.table(
      column = column_name,
      n_run_protein_groups = nrow(tmp),
      n_groups_with_multiple_values = sum(tmp$n_unique > 1L),
      max_unique_values = if (nrow(tmp) > 0L) max(tmp$n_unique) else NA_integer_
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

protein_output <- file.path(output_dir, "diann_protein_quant.tsv.gz")
validation_output <- file.path(output_dir, "diann_protein_field_consistency.tsv")
summary_output <- file.path(output_dir, "diann_extract_summary.txt")

data.table::fwrite(protein_quant, protein_output, sep = "\t", compress = "auto")
data.table::fwrite(consistency, validation_output, sep = "\t")

summary_lines <- c(
  paste0("Input file: ", input_file),
  sprintf("Input size (MiB): %.1f", file.info(input_file)$size / 1024^2),
  paste0("Selected columns: ", paste(available_columns, collapse = ", ")),
  paste0(
    "Protein quantity columns present: ",
    if (length(quant_columns) > 0L) paste(quant_columns, collapse = ", ") else "none"
  ),
  paste0("Report rows read: ", nrow(dt)),
  paste0("Runs: ", data.table::uniqueN(dt$Run)),
  paste0("Protein groups: ", data.table::uniqueN(dt$Protein.Group)),
  paste0("Run x protein groups: ", nrow(protein_quant)),
  "",
  "Outputs:",
  paste0("- ", file.path(output_dir, "diann_columns.txt")),
  paste0("- ", protein_output),
  paste0("- ", validation_output)
)
writeLines(summary_lines, summary_output)

if (nrow(consistency) > 0L && any(consistency$n_groups_with_multiple_values > 0L)) {
  warning(
    "At least one protein-level DIA-NN field has multiple values within a ",
    "Run x Protein.Group. Inspect diann_protein_field_consistency.tsv before ",
    "using the compact table."
  )
}

message("DIA-NN extraction complete.")
message("Summary: ", summary_output)
message("Protein table: ", protein_output)
