prepare_serum_metadata <- function(path, assay = c("massspec", "olink")) {
  assay <- match.arg(assay)

  raw <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  names(raw) <- trimws(names(raw))

  required <- c(
    "Unique sample name",
    "Exakt tube label (1)",
    "Group",
    "MassSpec",
    "Olink"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Raw metadata is missing columns: ", paste(missing, collapse = ", "))
  }

  assay_col <- if (assay == "massspec") "MassSpec" else "Olink"
  keep <- tolower(trimws(as.character(raw[[assay_col]]))) == "x"
  keep[is.na(keep)] <- FALSE
  dat <- raw[keep, , drop = FALSE]

  sample_id <- trimws(as.character(dat[["Unique sample name"]]))
  animal_id <- trimws(as.character(dat[["Exakt tube label (1)"]]))
  condition <- trimws(as.character(dat[["Group"]]))
  condition[condition == ""] <- NA_character_

  day_match <- regexec("_Day([0-9]+)$", sample_id)
  day_parts <- regmatches(sample_id, day_match)
  day <- vapply(
    day_parts,
    function(x) if (length(x) >= 2L) as.integer(x[[2L]]) else NA_integer_,
    integer(1)
  )

  sample_type <- ifelse(is.na(condition), "process_qc", "biological")
  visit <- ifelse(
    sample_type == "biological" & day == 0L,
    "baseline",
    ifelse(sample_type == "biological" & day > 0L, "endpoint", NA_character_)
  )

  out <- data.frame(
    sample_id = sample_id,
    animal_id = animal_id,
    condition = condition,
    day = day,
    visit = visit,
    sample_type = sample_type,
    stringsAsFactors = FALSE
  )

  endpoint_lookup <- out |>
    dplyr::filter(.data$sample_type == "biological") |>
    dplyr::group_by(.data$animal_id) |>
    dplyr::summarise(
      endpoint_day = {
        positive_days <- .data$day[!is.na(.data$day) & .data$day > 0L]
        if (length(positive_days) == 0L) NA_integer_ else max(positive_days)
      },
      .groups = "drop"
    )

  out <- dplyr::left_join(out, endpoint_lookup, by = "animal_id")
  out$condition_code <- gsub("[^A-Za-z0-9]+", "_", out$condition)
  out$condition_code <- gsub("^_|_$", "", out$condition_code)
  out$analysis_cell <- ifelse(
    out$sample_type == "biological",
    paste0(out$condition_code, "_D", out$endpoint_day, "_", out$visit),
    "ProcessQC"
  )

  if (anyDuplicated(out$sample_id)) {
    stop("Prepared metadata contains duplicated sample IDs.")
  }

  out |>
    dplyr::select(
      .data$sample_id,
      .data$animal_id,
      .data$condition,
      .data$condition_code,
      .data$day,
      .data$endpoint_day,
      .data$visit,
      .data$analysis_cell,
      .data$sample_type
    )
}

write_serum_metadata <- function(raw_metadata_path, output_dir = "metadata") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ms <- prepare_serum_metadata(raw_metadata_path, "massspec")
  olink <- prepare_serum_metadata(raw_metadata_path, "olink")

  readr::write_csv(ms, file.path(output_dir, "massspec_samples.csv"))
  readr::write_csv(olink, file.path(output_dir, "olink_samples.csv"))

  invisible(list(massspec = ms, olink = olink))
}
