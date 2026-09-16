`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

read_olink_config <- function(path = "config/olink.yml") {
  if (!file.exists(path)) stop("Olink config not found: ", path)
  yaml::read_yaml(path)
}

read_olink_paths <- function(path = "config/paths.yml") {
  if (!file.exists(path)) {
    stop("Copy config/paths.example.yml to config/paths.yml and edit local paths.")
  }
  paths <- yaml::read_yaml(path)
  if (is.null(paths$olink$quant) || is.null(paths$olink$metadata)) {
    stop("config/paths.yml must define olink.quant and olink.metadata.")
  }
  paths$olink
}

find_olink_row <- function(first_column, label) {
  idx <- which(trimws(as.character(first_column)) == label)
  if (length(idx) == 0L) stop("Could not find Olink row labelled '", label, "'.")
  idx[[1L]]
}

read_olink_quant <- function(path, metadata_path, cfg) {
  raw <- readxl::read_excel(
    path,
    sheet = cfg$input$sheet %||% 1,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  raw <- as.data.frame(raw, check.names = FALSE)
  names(raw) <- paste0("V", seq_len(ncol(raw)))

  first <- raw[[1L]]
  assay_row <- find_olink_row(first, "Assay")
  uniprot_row <- find_olink_row(first, "Uniprot ID")
  olinkid_row <- find_olink_row(first, "OlinkID")
  unit_row <- find_olink_row(first, "Unit")

  assay_values <- trimws(as.character(raw[assay_row, , drop = TRUE]))
  plate_col <- which(assay_values == "Plate ID")
  qc_col <- which(assay_values == "QC Warning")
  assay_cols <- which(
    seq_along(assay_values) > 1L &
      !is.na(assay_values) &
      nzchar(assay_values) &
      !assay_values %in% c("Plate ID", "QC Warning")
  )

  if (length(assay_cols) == 0L) stop("No Olink assay columns found.")
  if (length(qc_col) != 1L) stop("Expected one 'QC Warning' column in Olink workbook.")

  annotations <- data.frame(
    assay = assay_values[assay_cols],
    uniprot_id = as.character(raw[uniprot_row, assay_cols, drop = TRUE]),
    olink_id = as.character(raw[olinkid_row, assay_cols, drop = TRUE]),
    unit = as.character(raw[unit_row, assay_cols, drop = TRUE]),
    stringsAsFactors = FALSE
  )

  candidate_rows <- seq.int(unit_row + 1L, nrow(raw))
  sample_ids <- trimws(as.character(raw[candidate_rows, 1L, drop = TRUE]))
  keep <- !is.na(sample_ids) & nzchar(sample_ids)
  data_rows <- candidate_rows[keep]
  sample_ids <- sample_ids[keep]

  concentration <- raw[data_rows, assay_cols, drop = FALSE]
  concentration[] <- lapply(
    concentration,
    function(x) suppressWarnings(as.numeric(as.character(x)))
  )
  names(concentration) <- annotations$assay

  wide <- data.frame(sample_id = sample_ids, concentration, check.names = FALSE)
  wide$plate_id <- if (length(plate_col) == 1L) {
    as.character(raw[data_rows, plate_col, drop = TRUE])
  } else {
    NA_character_
  }
  wide$qc_warning <- as.character(raw[data_rows, qc_col, drop = TRUE])

  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)
  metadata$sample_id <- trimws(as.character(metadata$sample_id))

  unknown <- setdiff(wide$sample_id, metadata$sample_id)
  if (length(unknown) > 0L) {
    stop("Olink workbook samples absent from metadata: ", paste(unknown, collapse = ", "))
  }
  absent <- setdiff(metadata$sample_id, wide$sample_id)
  if (length(absent) > 0L) {
    warning("Olink metadata samples absent from workbook: ", paste(absent, collapse = ", "))
  }

  long <- wide |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(annotations$assay),
      names_to = "assay",
      values_to = "concentration"
    ) |>
    dplyr::left_join(annotations, by = "assay") |>
    dplyr::left_join(metadata, by = "sample_id")

  list(long = long, annotations = annotations, metadata = metadata)
}

prepare_olink_data <- function(long, cfg) {
  out <- long

  if (isTRUE(cfg$preprocessing$exclude_failed_qc)) {
    out <- out |>
      dplyr::filter(is.na(.data$qc_warning) | .data$qc_warning == "Pass")
  }

  if (isTRUE(cfg$preprocessing$zero_or_negative_to_na)) {
    out$concentration[!is.na(out$concentration) & out$concentration <= 0] <- NA_real_
  }

  if (isTRUE(cfg$preprocessing$impute_missing)) {
    stop("Olink missing-value imputation is intentionally not implemented in this scaffold.")
  }

  if (isTRUE(cfg$preprocessing$log2_transform)) {
    out$log2_concentration <- log2(out$concentration)
  } else {
    out$log2_concentration <- out$concentration
  }

  out$animal_id <- factor(out$animal_id)
  out$analysis_cell <- factor(out$analysis_cell)
  out$sample_type <- factor(out$sample_type)
  out
}

olink_qc_tables <- function(long) {
  sample_qc <- long |>
    dplyr::group_by(.data$sample_id, .data$qc_warning, .data$plate_id) |>
    dplyr::summarise(
      n_assays = dplyr::n(),
      n_quantified = sum(!is.na(.data$concentration)),
      missing_fraction = mean(is.na(.data$concentration)),
      median_pg_ml = stats::median(.data$concentration, na.rm = TRUE),
      .groups = "drop"
    )

  assay_qc <- long |>
    dplyr::group_by(.data$assay, .data$uniprot_id, .data$olink_id, .data$unit) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_quantified = sum(!is.na(.data$concentration)),
      missing_fraction = mean(is.na(.data$concentration)),
      median_pg_ml = stats::median(.data$concentration, na.rm = TRUE),
      .groups = "drop"
    )

  design <- long |>
    dplyr::distinct(
      .data$sample_id,
      .data$animal_id,
      .data$condition,
      .data$endpoint_day,
      .data$visit,
      .data$analysis_cell,
      .data$sample_type
    ) |>
    dplyr::count(
      .data$condition,
      .data$endpoint_day,
      .data$visit,
      .data$sample_type,
      name = "n_samples"
    )

  list(sample = sample_qc, assay = assay_qc, design = design)
}

olink_design_columns <- function(long, cfg) {
  dat <- long |>
    dplyr::filter(.data$sample_type == "biological") |>
    dplyr::distinct(.data$sample_id, .keep_all = TRUE)
  fixed <- cfg$model$fixed_effects %||% "0 + analysis_cell"
  colnames(stats::model.matrix(stats::as.formula(paste("~", fixed)), data = dat))
}

fit_olink_contrasts <- function(long, cfg) {
  contrasts <- unlist(cfg$model$contrasts, use.names = TRUE)
  if (length(contrasts) == 0L) {
    stop(
      "No Olink contrasts configured. Inspect results/proteomics/olink/model_matrix_columns.txt ",
      "and add named expressions under model.contrasts in config/olink.yml."
    )
  }

  dat <- long |>
    dplyr::filter(.data$sample_type == "biological")

  sample_meta <- dat |>
    dplyr::distinct(.data$sample_id, .keep_all = TRUE) |>
    dplyr::arrange(.data$sample_id)

  fixed <- cfg$model$fixed_effects %||% "0 + analysis_cell"
  design <- stats::model.matrix(stats::as.formula(paste("~", fixed)), data = sample_meta)
  if (qr(design)$rank < ncol(design)) {
    stop("Olink model matrix is rank-deficient. Inspect the observed design cells before defining contrasts.")
  }

  duplicate_pairs <- dat |>
    dplyr::count(.data$assay, .data$sample_id) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicate_pairs) > 0L) {
    stop("Olink data contain duplicated assay x sample measurements.")
  }

  wide <- dat |>
    dplyr::select(.data$assay, .data$sample_id, .data$log2_concentration) |>
    tidyr::pivot_wider(names_from = .data$sample_id, values_from = .data$log2_concentration)

  sample_ids <- sample_meta$sample_id
  missing_samples <- setdiff(sample_ids, names(wide))
  if (length(missing_samples) > 0L) {
    stop("Olink matrix is missing samples: ", paste(missing_samples, collapse = ", "))
  }

  mat <- as.matrix(wide[, sample_ids, drop = FALSE])
  rownames(mat) <- wide$assay
  storage.mode(mat) <- "double"

  max_missing <- cfg$filtering$max_missing_fraction
  if (!is.null(max_missing)) {
    keep <- rowMeans(is.na(mat)) <= as.numeric(max_missing)
    mat <- mat[keep, , drop = FALSE]
  }

  contrast_matrix <- limma::makeContrasts(
    contrasts = unname(contrasts),
    levels = design
  )
  colnames(contrast_matrix) <- names(contrasts)

  animal <- sample_meta$animal_id
  repeated <- any(table(animal) > 1L)
  correlation <- NA_real_

  if (repeated && identical(cfg$model$repeated_measure_method, "duplicateCorrelation")) {
    corfit <- limma::duplicateCorrelation(mat, design, block = animal)
    correlation <- corfit$consensus.correlation
    if (!is.finite(correlation)) {
      stop("limma duplicateCorrelation did not return a finite consensus correlation.")
    }
    fit <- limma::lmFit(mat, design, block = animal, correlation = correlation)
  } else {
    fit <- limma::lmFit(mat, design)
  }

  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(
    fit,
    trend = isTRUE(cfg$model$limma$trend),
    robust = isTRUE(cfg$model$limma$robust)
  )

  results <- lapply(seq_len(ncol(contrast_matrix)), function(i) {
    tab <- limma::topTable(fit, coef = i, number = Inf, sort.by = "none")
    tab$assay <- rownames(tab)
    tab$contrast <- colnames(contrast_matrix)[[i]]
    tibble::as_tibble(tab)
  }) |>
    dplyr::bind_rows() |>
    dplyr::select(.data$assay, .data$contrast, dplyr::everything())

  list(
    design = design,
    contrast_matrix = contrast_matrix,
    correlation = correlation,
    fit = fit,
    results = results
  )
}
