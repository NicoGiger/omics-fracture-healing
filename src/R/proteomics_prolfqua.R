`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

read_proteomics_config <- function(path = "config/proteomics.yml") {
  if (!file.exists(path)) {
    stop("Proteomics config not found: ", path)
  }
  yaml::read_yaml(path)
}

read_proteomics_paths <- function(path = "config/paths.yml") {
  if (!file.exists(path)) {
    stop(
      "Local path config not found: ", path,
      ". Copy config/paths.example.yml to config/paths.yml and edit it."
    )
  }

  paths <- yaml::read_yaml(path)
  if (is.null(paths$proteomics$abundance) || is.null(paths$proteomics$metadata)) {
    stop("config/paths.yml must define proteomics.abundance and proteomics.metadata.")
  }
  paths$proteomics
}

validate_proteomics_metadata <- function(metadata, cfg) {
  cols <- cfg$columns
  required <- c(cols$sample_id, cols$animal_id, cols$condition, cols$timepoint)
  required <- required[!vapply(required, is.null, logical(1))]

  missing_cols <- setdiff(required, names(metadata))
  if (length(missing_cols) > 0L) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  sample_id <- cols$sample_id
  if (anyNA(metadata[[sample_id]]) || any(metadata[[sample_id]] == "")) {
    stop("Metadata contains missing/empty sample IDs.")
  }
  if (anyDuplicated(metadata[[sample_id]])) {
    stop("Metadata sample IDs must be unique.")
  }

  invisible(metadata)
}

coerce_proteomics_factors <- function(data, cfg) {
  cols <- cfg$columns

  if (!is.null(cols$condition) && isTRUE(cfg$model$condition_as_factor %||% TRUE)) {
    data[[cols$condition]] <- factor(data[[cols$condition]])
  }
  if (!is.null(cols$timepoint) && isTRUE(cfg$model$timepoint_as_factor %||% TRUE)) {
    data[[cols$timepoint]] <- factor(data[[cols$timepoint]])
  }
  if (!is.null(cols$animal_id)) {
    data[[cols$animal_id]] <- factor(data[[cols$animal_id]])
  }
  if (!is.null(cols$batch) && isTRUE(cfg$model$batch_as_factor %||% TRUE)) {
    data[[cols$batch]] <- factor(data[[cols$batch]])
  }

  data
}

read_proteomics_wide <- function(abundance_path, metadata_path, cfg) {
  if (!identical(cfg$input$format, "wide_protein")) {
    stop("Current adapter supports input.format = 'wide_protein' only.")
  }

  delimiter <- cfg$input$delimiter %||% "\t"
  abundance <- readr::read_delim(
    abundance_path,
    delim = delimiter,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)

  validate_proteomics_metadata(metadata, cfg)

  protein_id <- cfg$columns$protein_id
  sample_id <- cfg$columns$sample_id

  if (!protein_id %in% names(abundance)) {
    stop("Protein ID column not found in abundance table: ", protein_id)
  }

  sample_ids <- as.character(metadata[[sample_id]])
  absent_samples <- setdiff(sample_ids, names(abundance))
  if (length(absent_samples) > 0L) {
    stop(
      "Samples present in metadata but absent from abundance table: ",
      paste(absent_samples, collapse = ", ")
    )
  }

  long <- abundance |>
    dplyr::select(dplyr::all_of(c(protein_id, sample_ids))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_ids),
      names_to = sample_id,
      values_to = "intensity"
    ) |>
    dplyr::left_join(metadata, by = sample_id)

  long[[protein_id]] <- as.character(long[[protein_id]])
  long[[sample_id]] <- as.character(long[[sample_id]])
  long <- coerce_proteomics_factors(long, cfg)

  if (isTRUE(cfg$preprocessing$zero_or_negative_to_na)) {
    long$intensity[!is.na(long$intensity) & long$intensity <= 0] <- NA_real_
  }

  list(data = long, metadata = metadata)
}

make_prolfqua_data <- function(long_data, cfg) {
  cols <- cfg$columns

  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- cols$sample_id
  config$work_intensity <- "intensity"
  config$hierarchy[[cols$protein_id]] <- cols$protein_id

  factor_cols <- unique(c(
    cols$condition,
    cols$timepoint,
    cols$animal_id,
    cols$batch
  ))
  factor_cols <- factor_cols[
    !is.na(factor_cols) & nzchar(factor_cols)
  ]

  for (factor_col in factor_cols) {
    if (!factor_col %in% names(long_data)) {
      stop("Configured factor column is absent from data: ", factor_col)
    }
    config$factors[[factor_col]] <- factor_col
  }

  analysis_data <- prolfqua::setup_analysis(long_data, config)
  prolfqua::LFQData$new(analysis_data, config)
}

preprocess_prolfqua <- function(lfq, cfg) {
  scale <- cfg$preprocessing$intensity_scale %||% "linear"
  normalization <- cfg$preprocessing$normalization %||% "none"

  if (!scale %in% c("linear", "log2")) {
    stop("preprocessing.intensity_scale must be 'linear' or 'log2'.")
  }
  if (!normalization %in% c("none", "robscale")) {
    stop("preprocessing.normalization must be 'none' or 'robscale'.")
  }

  transformer <- lfq$get_Transformer()
  if (identical(scale, "linear")) {
    transformer <- transformer$log2()
  }
  if (identical(normalization, "robscale")) {
    transformer <- transformer$robscale()
  }

  out <- transformer$lfq

  response_name <- if (identical(normalization, "robscale")) {
    "normalized_log2_intensity"
  } else if (identical(scale, "linear")) {
    "log2_intensity"
  } else {
    "intensity"
  }

  out$rename_response(response_name)
  out
}

proteomics_qc_tables <- function(lfq, cfg) {
  dat <- lfq$data_long()
  response <- lfq$response()
  sample_id <- cfg$columns$sample_id
  protein_id <- cfg$columns$protein_id

  sample_qc <- dat |>
    dplyr::group_by(.data[[sample_id]]) |>
    dplyr::summarise(
      n_proteins = dplyr::n(),
      n_detected = sum(!is.na(.data[[response]])),
      missing_fraction = mean(is.na(.data[[response]])),
      median_intensity = stats::median(.data[[response]], na.rm = TRUE),
      .groups = "drop"
    )

  protein_qc <- dat |>
    dplyr::group_by(.data[[protein_id]]) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_detected = sum(!is.na(.data[[response]])),
      missing_fraction = mean(is.na(.data[[response]])),
      median_intensity = stats::median(.data[[response]], na.rm = TRUE),
      .groups = "drop"
    )

  list(sample = sample_qc, protein = protein_qc)
}

has_repeated_animals <- function(lfq, cfg) {
  dat <- lfq$data_long()
  sample_id <- cfg$columns$sample_id
  animal_id <- cfg$columns$animal_id

  sample_map <- dat |>
    dplyr::distinct(.data[[sample_id]], .data[[animal_id]])

  any(table(sample_map[[animal_id]]) > 1L)
}

proteomics_fixed_effects <- function(cfg) {
  fixed <- cfg$model$fixed_effects %||% "auto"
  if (!identical(fixed, "auto")) {
    return(fixed)
  }

  condition <- cfg$columns$condition
  timepoint <- cfg$columns$timepoint
  batch <- cfg$columns$batch

  fixed <- if (!is.null(condition) && !is.null(timepoint)) {
    paste(condition, "*", timepoint)
  } else if (!is.null(condition)) {
    condition
  } else if (!is.null(timepoint)) {
    timepoint
  } else {
    stop("No fixed effects are configured.")
  }

  if (!is.null(batch)) {
    fixed <- paste(fixed, "+", batch)
  }

  fixed
}

proteomics_design_columns <- function(lfq, cfg) {
  dat <- lfq$data_long()
  sample_id <- cfg$columns$sample_id
  fixed <- proteomics_fixed_effects(cfg)

  sample_data <- dat |>
    dplyr::distinct(.data[[sample_id]], .keep_all = TRUE)

  mm <- stats::model.matrix(stats::as.formula(paste("~", fixed)), data = sample_data)
  colnames(mm)
}

fit_proteomics_contrasts <- function(lfq, cfg) {
  contrasts <- unlist(cfg$model$contrasts, use.names = TRUE)
  if (length(contrasts) == 0L) {
    stop(
      "No contrasts configured. Inspect results/proteomics/model_matrix_columns.txt ",
      "and then add named expressions under model.contrasts in config/proteomics.yml."
    )
  }

  mode <- cfg$model$mode %||% "auto"
  repeated <- has_repeated_animals(lfq, cfg)

  if (identical(mode, "auto")) {
    mode <- if (repeated) "lmer" else "limma"
  }

  if (!mode %in% c("limma", "lmer")) {
    stop("model.mode must be one of: auto, limma, lmer.")
  }

  if (identical(mode, "limma") && repeated) {
    warning(
      "Repeated animals were detected but model.mode = 'limma'. ",
      "This treats samples as independent. Use only if this is intentional."
    )
  }

  response <- lfq$response()
  fixed <- proteomics_fixed_effects(cfg)

  if (identical(mode, "limma")) {
    formula <- paste(response, "~", fixed)
    strat <- prolfqua::strategy_limma(
      formula,
      trend = isTRUE(cfg$model$limma$trend),
      robust = isTRUE(cfg$model$limma$robust)
    )
    model <- prolfqua::build_model_limma(lfq, strat)
    contrast_object <- prolfqua::ContrastsLimma$new(model, contrasts)
  } else {
    animal_id <- cfg$columns$animal_id
    formula <- paste(response, "~", fixed, "+ (1 |", animal_id, ")")
    strat <- prolfqua::strategy_lmer(formula)
    model <- prolfqua::build_model(lfq, strat)
    contrast_object <- prolfqua::Contrasts$new(model, contrasts)
  }

  results <- contrast_object$get_contrasts()
  if ("p.value" %in% names(results) && "contrast" %in% names(results)) {
    results <- results |>
      dplyr::group_by(.data$contrast) |>
      dplyr::mutate(FDR_BH = stats::p.adjust(.data$p.value, method = "BH")) |>
      dplyr::ungroup()
  }

  list(
    mode = mode,
    formula = formula,
    model = model,
    contrasts = contrast_object,
    results = results
  )
}

save_prolfqua_bundle <- function(lfq, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  bundle <- list(
    data = lfq$data_long(),
    config = prolfqua::R6_extract_values(lfq$get_config())
  )
  saveRDS(bundle, path)
  invisible(path)
}

load_prolfqua_bundle <- function(path) {
  bundle <- readRDS(path)
  config <- prolfqua::list_to_AnalysisConfiguration(bundle$config)
  prolfqua::LFQData$new(bundle$data, config)
}
