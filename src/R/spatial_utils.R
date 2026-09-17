read_spatial_samples <- function(path = "metadata/spatial_samples.csv") {
  samples <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )

  required <- c(
    "capture_id", "time_point", "bio_rep", "tech_rep", "data_subdir"
  )
  missing <- setdiff(required, names(samples))
  if (length(missing) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(samples$capture_id)) {
    stop("capture_id must be unique in ", path)
  }

  samples$bio_rep <- as.integer(samples$bio_rep)
  samples$tech_rep <- as.integer(samples$tech_rep)
  samples
}

read_spot_annotations <- function(path) {
  if (!file.exists(path)) stop("Annotation file not found: ", path)

  annotation <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  names_lower <- tolower(names(annotation))
  barcode_idx <- match("barcode", names_lower)
  if (is.na(barcode_idx)) stop("Expected a Barcode column in ", path)

  barcodes <- annotation[[barcode_idx]]
  if (anyDuplicated(barcodes)) {
    stop("Duplicate barcodes in annotation file: ", path)
  }

  get_col <- function(candidates) {
    idx <- match(tolower(candidates), names_lower, nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) {
      rep(NA_character_, nrow(annotation))
    } else {
      annotation[[idx[[1L]]]]
    }
  }

  data.frame(
    condition = get_col(c("condition", "group")),
    roi = get_col("roi"),
    tissue = get_col("tissue"),
    artifact = get_col(c("artifact", "artefact", "artifacts", "artefacts")),
    stringsAsFactors = FALSE,
    row.names = barcodes
  )
}

add_spot_annotations <- function(obj, path) {
  annotation <- read_spot_annotations(path)
  mapped <- annotation[match(colnames(obj), rownames(annotation)), , drop = FALSE]
  rownames(mapped) <- colnames(obj)
  SeuratObject::AddMetaData(obj, mapped)
}

artifact_is_flagged <- function(x) {
  if (all(is.na(x))) return(rep(FALSE, length(x)))
  x <- tolower(trimws(as.character(x)))
  !(is.na(x) | x == "" | x %in% c("no", "none", "false", "0"))
}

safe_median <- function(x) {
  if (length(x) == 0L) NA_real_ else stats::median(x, na.rm = TRUE)
}

summarize_condition_qc <- function(meta_before, meta_after, capture) {
  conditions <- sort(unique(meta_before$condition))

  do.call(rbind, lapply(conditions, function(condition) {
    before <- meta_before[meta_before$condition == condition, , drop = FALSE]
    after <- meta_after[meta_after$condition == condition, , drop = FALSE]

    data.frame(
      capture_id = capture$capture_id,
      sample_id = paste0(capture$time_point, "_", condition, "_r", capture$bio_rep),
      section_id = paste(capture$capture_id, condition, sep = "_"),
      time_point = capture$time_point,
      condition = condition,
      bio_rep = capture$bio_rep,
      tech_rep = capture$tech_rep,
      n_spots_before_qc = nrow(before),
      n_spots_retained = nrow(after),
      retained_fraction = if (nrow(before) > 0L) nrow(after) / nrow(before) else NA_real_,
      median_total_counts = safe_median(after$nCount_Spatial),
      median_detected_genes = safe_median(after$nFeature_Spatial),
      median_percent_mt = safe_median(after$percent.mt),
      row.names = NULL
    )
  }))
}

process_spatial_capture <- function(
    capture,
    data_root,
    min_features = 100L,
    max_percent_mt = 10
) {
  data_dir <- file.path(data_root, capture$data_subdir)
  counts_file <- file.path(data_dir, "raw_feature_bc_matrix.h5")
  annotation_file <- file.path(data_dir, "spa.csv")

  if (!dir.exists(data_dir)) stop("Spatial data directory not found: ", data_dir)
  if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
  if (!file.exists(annotation_file)) stop("Annotation file not found: ", annotation_file)

  obj <- Seurat::Load10X_Spatial(
    data.dir = data_dir,
    filename = basename(counts_file),
    assay = "Spatial",
    filter.matrix = FALSE
  )
  n_loaded <- ncol(obj)

  obj <- add_spot_annotations(obj, annotation_file)

  condition <- trimws(as.character(obj$condition))
  condition[tolower(condition) == "union"] <- "Union"
  condition[tolower(condition) == "nonunion"] <- "NonUnion"
  obj$condition <- condition

  condition_assigned <- !is.na(obj$condition) & obj$condition %in% c("Union", "NonUnion")
  n_condition_assigned <- sum(condition_assigned)
  if (n_condition_assigned == 0L) {
    stop("No Union/NonUnion spots assigned for ", capture$capture_id)
  }

  n_roi_assigned <- sum(condition_assigned & !is.na(obj$roi) & obj$roi != "")
  n_tissue_assigned <- sum(condition_assigned & !is.na(obj$tissue) & obj$tissue != "")
  tissue_available <- n_tissue_assigned > 0L

  obj <- subset(obj, cells = colnames(obj)[condition_assigned])

  artifact_flag <- artifact_is_flagged(obj$artifact)
  n_artifact_flagged <- sum(artifact_flag)
  if (n_artifact_flagged > 0L) {
    obj <- subset(obj, cells = colnames(obj)[!artifact_flag])
  }
  if (ncol(obj) == 0L) {
    stop("No spots left after annotation/artifact filtering for ", capture$capture_id)
  }

  mt_genes <- grep("^mt-", rownames(obj), ignore.case = TRUE, value = TRUE)
  if (length(mt_genes) == 0L) {
    warning("No mitochondrial genes detected in ", capture$capture_id)
    obj$percent.mt <- 0
  } else {
    obj <- Seurat::PercentageFeatureSet(
      obj,
      features = mt_genes,
      col.name = "percent.mt"
    )
  }

  obj$capture_id <- capture$capture_id
  obj$time_point <- capture$time_point
  obj$bio_rep <- capture$bio_rep
  obj$tech_rep <- capture$tech_rep
  obj$sample_id <- paste0(obj$time_point, "_", obj$condition, "_r", obj$bio_rep)
  obj$section_id <- paste(obj$capture_id, obj$condition, sep = "_")
  obj$orig.ident <- obj$section_id

  meta_before_qc <- obj[[]]
  keep <- meta_before_qc$nFeature_Spatial >= min_features &
    meta_before_qc$percent.mt < max_percent_mt
  keep[is.na(keep)] <- FALSE
  obj <- subset(obj, cells = rownames(meta_before_qc)[keep])
  if (ncol(obj) == 0L) {
    stop("No spots left after QC for ", capture$capture_id)
  }
  meta_after_qc <- obj[[]]

  annotation_coverage <- data.frame(
    capture_id = capture$capture_id,
    annotation_file = basename(annotation_file),
    time_point = capture$time_point,
    bio_rep = capture$bio_rep,
    tech_rep = capture$tech_rep,
    n_barcodes_loaded = n_loaded,
    n_condition_assigned = n_condition_assigned,
    n_roi_assigned = n_roi_assigned,
    roi_assignment_fraction = n_roi_assigned / n_condition_assigned,
    tissue_annotation_available = tissue_available,
    n_tissue_assigned = n_tissue_assigned,
    tissue_assignment_fraction = n_tissue_assigned / n_condition_assigned,
    n_artifact_flagged = n_artifact_flagged,
    n_spots_before_qc = nrow(meta_before_qc),
    n_spots_retained = nrow(meta_after_qc),
    row.names = NULL
  )

  list(
    object = obj,
    sample_qc = summarize_condition_qc(meta_before_qc, meta_after_qc, capture),
    annotation_coverage = annotation_coverage
  )
}

save_spatial_qc_plots <- function(obj, capture_id, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(
    file.path(output_dir, paste0(capture_id, "_qc.pdf")),
    width = 8,
    height = 7
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  print(
    Seurat::SpatialDimPlot(obj, group.by = "condition") +
      ggplot2::ggtitle(paste(capture_id, "condition"))
  )
  print(
    Seurat::SpatialFeaturePlot(obj, features = "nCount_Spatial") +
      ggplot2::ggtitle(paste(capture_id, "total counts"))
  )
  print(
    Seurat::SpatialFeaturePlot(obj, features = "nFeature_Spatial") +
      ggplot2::ggtitle(paste(capture_id, "detected genes"))
  )
  print(
    Seurat::SpatialFeaturePlot(obj, features = "percent.mt") +
      ggplot2::ggtitle(paste(capture_id, "mitochondrial fraction"))
  )

  if (any(!is.na(obj$roi) & obj$roi != "")) {
    print(
      Seurat::SpatialDimPlot(obj, group.by = "roi") +
        ggplot2::ggtitle(paste(capture_id, "ROI"))
    )
  }
  if (any(!is.na(obj$tissue) & obj$tissue != "")) {
    print(
      Seurat::SpatialDimPlot(obj, group.by = "tissue") +
        ggplot2::ggtitle(paste(capture_id, "tissue"))
    )
  }

  invisible(NULL)
}
