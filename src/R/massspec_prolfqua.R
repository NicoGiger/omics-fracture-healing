`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

read_ms_config <- function(path = "config/massspec.yml") {
  if (!file.exists(path)) stop("Mass-spec config not found: ", path)
  yaml::read_yaml(path)
}

read_ms_paths <- function(path = "config/paths.yml") {
  if (!file.exists(path)) {
    stop("Copy config/paths.example.yml to config/paths.yml and edit local paths.")
  }
  paths <- yaml::read_yaml(path)
  if (is.null(paths$massspec$abundance) || is.null(paths$massspec$metadata)) {
    stop("config/paths.yml must define massspec.abundance and massspec.metadata.")
  }
  paths$massspec
}

read_ms_wide <- function(abundance_path, metadata_path, cfg) {
  delimiter <- cfg$input$delimiter %||% ","
  abundance <- readr::read_delim(
    abundance_path,
    delim = delimiter,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)

  cols <- cfg$columns
  required_meta <- c(cols$sample_id, cols$animal_id, cols$analysis_cell, cols$sample_type)
  missing_meta <- setdiff(required_meta, names(metadata))
  if (length(missing_meta) > 0L) {
    stop("Mass-spec metadata is missing columns: ", paste(missing_meta, collapse = ", "))
  }
  if (!cols$protein_id %in% names(abundance)) {
    stop("Protein ID column not found: ", cols$protein_id)
  }

  metadata[[cols$sample_id]] <- trimws(as.character(metadata[[cols$sample_id]]))
  abundance_samples <- setdiff(names(abundance), cols$protein_id)
  metadata_samples <- metadata[[cols$sample_id]]

  unknown_abundance <- setdiff(abundance_samples, metadata_samples)
  if (length(unknown_abundance) > 0L) {
    stop(
      "Abundance matrix contains samples absent from metadata: ",
      paste(unknown_abundance, collapse = ", ")
    )
  }

  missing_abundance <- setdiff(metadata_samples, abundance_samples)
  if (length(missing_abundance) > 0L) {
    warning(
      "Metadata samples absent from abundance matrix (kept in metadata only): ",
      paste(missing_abundance, collapse = ", ")
    )
  }

  long <- abundance |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(abundance_samples),
      names_to = cols$sample_id,
      values_to = "intensity"
    ) |>
    dplyr::left_join(metadata, by = cols$sample_id)

  long[[cols$protein_id]] <- as.character(long[[cols$protein_id]])
  long[[cols$sample_id]] <- as.character(long[[cols$sample_id]])
  long[[cols$animal_id]] <- factor(long[[cols$animal_id]])
  long[[cols$analysis_cell]] <- factor(long[[cols$analysis_cell]])
  long[[cols$sample_type]] <- factor(long[[cols$sample_type]])
  if (!is.null(cols$batch) && cols$batch %in% names(long)) {
    long[[cols$batch]] <- factor(long[[cols$batch]])
  }

  if (isTRUE(cfg$preprocessing$zero_or_negative_to_na)) {
    long$intensity[!is.na(long$intensity) & long$intensity <= 0] <- NA_real_
  }

  list(data = long, metadata = metadata)
}

make_ms_lfq <- function(long_data, cfg) {
  cols <- cfg$columns

  # Standard prolfqua entry point: AnalysisConfiguration -> setup_analysis -> LFQData.
  conf <- prolfqua::AnalysisConfiguration$new()
  conf$file_name <- cols$sample_id
  conf$work_intensity <- "intensity"
  conf$hierarchy[[cols$protein_id]] <- cols$protein_id

  factor_cols <- unique(c(
    cols$animal_id,
    cols$analysis_cell,
    cols$sample_type,
    cols$batch
  ))
  factor_cols <- factor_cols[!is.na(factor_cols) & nzchar(factor_cols)]
  for (factor_col in factor_cols) {
    conf$factors[[factor_col]] <- factor_col
  }

  analysed <- prolfqua::setup_analysis(long_data, conf)
  prolfqua::LFQData$new(analysed, conf)
}

preprocess_ms_lfq <- function(lfq, cfg) {
  scale <- cfg$preprocessing$intensity_scale %||% "linear"
  normalization <- cfg$preprocessing$normalization %||% "none"

  if (!scale %in% c("linear", "log2")) {
    stop("massspec preprocessing.intensity_scale must be linear or log2.")
  }
  if (!normalization %in% c("none", "robscale")) {
    stop("massspec preprocessing.normalization must be none or robscale.")
  }

  transformer <- lfq$get_Transformer()
  if (scale == "linear") transformer <- transformer$log2()
  if (normalization == "robscale") transformer <- transformer$robscale()

  out <- transformer$lfq
  response_name <- if (normalization == "robscale") {
    "normalized_log2_intensity"
  } else if (scale == "linear") {
    "log2_intensity"
  } else {
    "intensity"
  }
  out$rename_response(response_name)
  out
}

ms_qc_tables <- function(lfq, cfg) {
  dat <- lfq$data_long()
  response <- lfq$response()
  sample_id <- cfg$columns$sample_id
  protein_id <- cfg$columns$protein_id
  sample_type <- cfg$columns$sample_type

  sample_qc <- dat |>
    dplyr::group_by(.data[[sample_id]], .data[[sample_type]]) |>
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

ms_biological_lfq <- function(lfq, cfg) {
  sample_type <- cfg$columns$sample_type
  dat <- lfq$data_long() |>
    dplyr::filter(.data[[sample_type]] == "biological")
  prolfqua::LFQData$new(dat, lfq$get_config())
}

ms_design_columns <- function(lfq, cfg) {
  cols <- cfg$columns
  dat <- lfq$data_long() |>
    dplyr::distinct(.data[[cols$sample_id]], .keep_all = TRUE)
  fixed <- cfg$model$fixed_effects %||% "0 + analysis_cell"
  mm <- stats::model.matrix(stats::as.formula(paste("~", fixed)), data = dat)
  colnames(mm)
}

has_repeated_ms_animals <- function(lfq, cfg) {
  cols <- cfg$columns
  dat <- lfq$data_long() |>
    dplyr::distinct(.data[[cols$sample_id]], .data[[cols$animal_id]])
  any(table(dat[[cols$animal_id]]) > 1L)
}

fit_ms_contrasts <- function(lfq, cfg) {
  contrasts <- unlist(cfg$model$contrasts, use.names = TRUE)
  if (length(contrasts) == 0L) {
    stop(
      "No MS contrasts configured. Inspect results/proteomics/massspec/model_matrix_columns.txt ",
      "and add named expressions under model.contrasts in config/massspec.yml."
    )
  }

  mode <- cfg$model$mode %||% "auto"
  repeated <- has_repeated_ms_animals(lfq, cfg)
  if (mode == "auto") mode <- if (repeated) "lmer" else "limma"
  if (!mode %in% c("limma", "lmer")) stop("MS model.mode must be auto, limma, or lmer.")

  response <- lfq$response()
  fixed <- cfg$model$fixed_effects %||% "0 + analysis_cell"

  if (mode == "limma") {
    if (repeated) {
      warning("Repeated animals detected: plain prolfqua limma treats runs as independent.")
    }
    formula <- paste(response, "~", fixed)
    strategy <- prolfqua::strategy_limma(
      formula,
      trend = isTRUE(cfg$model$limma$trend),
      robust = isTRUE(cfg$model$limma$robust)
    )
    model <- prolfqua::build_model_limma(lfq, strategy)
    contrast_object <- prolfqua::ContrastsLimma$new(model, contrasts)
  } else {
    animal_id <- cfg$columns$animal_id
    formula <- paste(response, "~", fixed, "+ (1 |", animal_id, ")")
    strategy <- prolfqua::strategy_lmer(formula)
    model <- prolfqua::build_model(lfq, strategy)
    contrast_object <- prolfqua::Contrasts$new(model, contrasts)
  }

  results <- contrast_object$get_contrasts()
  if (all(c("contrast", "p.value") %in% names(results))) {
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

save_ms_lfq <- function(lfq, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  bundle <- list(
    data = lfq$data_long(),
    config = prolfqua::R6_extract_values(lfq$get_config())
  )
  saveRDS(bundle, path)
  invisible(path)
}

load_ms_lfq <- function(path) {
  bundle <- readRDS(path)
  conf <- prolfqua::list_to_AnalysisConfiguration(bundle$config)
  prolfqua::LFQData$new(bundle$data, conf)
}
