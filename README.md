# omics-fracture-healing

Paper-facing analysis repository for the fracture-healing study, integrating spatial transcriptomics and serum proteomics.

## Repository structure

- `scripts/` — executable analysis entry points. Keep these thin; reusable logic belongs in `src/`.
- `src/R/` — reusable R functions for metadata, spatial transcriptomics, serum proteomics, statistics, enrichment, integration, and plotting.
- `src/python/` — reusable Python code, including the SpatialDM pipeline.
- `metadata/` — distilled assay-specific metadata derived from the full study metadata. Spatial and serum tables should share a stable `animal_id` key.
- `config/` — analysis configuration. Machine-specific paths should live in `config/paths.yml` and remain untracked.
- `data/` — local/raw or intermediate data locations. Large source data should not be committed unless explicitly intended.
- `analysis/` — Quarto/R Markdown reports and paper-facing exploratory summaries that consume processed results.
- `results/` — derived numerical outputs from the analysis pipelines.
- `figures/` — exploratory, manuscript, and supplementary figures.

## Planned analysis branches

### Spatial transcriptomics

1. Metadata validation and raw-data loading
2. Spot-level QC with SpotSweeper
3. Cell-type deconvolution
4. Spot-level pathway activity
5. Spatial organization analyses (including Moran's I)
6. SpatialDM ligand-receptor analysis
7. Raw-count pseudobulk analysis
8. Temporal trajectory / shape-shift analysis
9. ROI-aware analyses and pathway enrichment

### Serum proteomics

The serum workflow uses [`prolfqua`](https://github.com/fgcz/prolfqua) for proteomics data representation, preprocessing/QC, and differential modelling. The environment setup pins the stable `prolfqua` v1.5.0 release in `renv` rather than tracking the moving development branch.

#### DIA-NN source inspection

The original DIA-NN report is large and should remain local. `scripts/00_extract_diann_quant.R` reads only the columns needed to determine which protein-level quantity and filtering rules should feed the analysis.

After initializing the R environment, run:

```bash
Rscript scripts/00_extract_diann_quant.R /path/to/report.tsv data/proteomics/diann_extract
```

The script does not alter or filter the source report. It writes:

- `diann_columns.txt` — all DIA-NN report column names.
- `diann_protein_quant.tsv.gz` — compact `Run x Protein.Group` table containing available protein-level quantities (`PG.MaxLFQ`, `PG.Normalised`, and/or `PG.Quantity`), protein-group Q-values, annotation fields, and precursor/Q-value counts.
- `diann_protein_field_consistency.tsv` — checks whether protein-level DIA-NN fields are actually constant across precursor rows for each `Run x Protein.Group`.
- `diann_extract_summary.txt` — source size, selected fields, run/protein counts, and output locations.

This compact extract is intended for local inspection before deciding whether the existing `PG.Normalised` wide matrix should remain the canonical mass-spectrometry input or be regenerated from another DIA-NN protein quantity. The raw DIA-NN report itself should not be committed.

#### Current prolfqua workflow

1. Copy `config/paths.example.yml` to `config/paths.yml` and set local abundance/metadata paths.
2. Create `metadata/proteomics_samples.csv` from `metadata/proteomics_samples_template.csv`.
3. Edit `config/proteomics.yml` to match column names, intensity scale, normalization choice, and the experimental design.
4. Run `Rscript scripts/00_setup_R.R` once to initialize/snapshot the R environment.
5. Optionally run the DIA-NN source inspection above before fixing the canonical mass-spec abundance input.
6. Run `Rscript scripts/01_proteomics_qc.R` to construct the `LFQData`, preprocess it, write QC tables/plots, and export `results/proteomics/model_matrix_columns.txt`.
7. Define named model contrasts under `model.contrasts` in `config/proteomics.yml` after inspecting the real design matrix.
8. Run `Rscript scripts/02_proteomics_differential.R` for protein-level differential analysis.

Defaults are intentionally conservative: no missing-value imputation and no extra normalization are applied unless requested in the configuration. `model.mode: auto` uses limma for independent samples and switches to a random-intercept mixed model when repeated samples from the same `animal_id` are detected. This prevents accidental pseudoreplication, but the final fixed-effects structure and contrasts still need to be checked against the actual serum sampling design.

Planned downstream steps:

1. Assay-specific QC and preprocessing
2. Protein-level temporal and outcome analyses
3. Pathway-level interpretation

### Integration

Cross-omics analyses are performed only after assay-specific analyses, using `animal_id` as the shared experimental-unit key where measurements are matched.

## Reproducibility

R dependencies are managed with `renv`. Python dependencies for SpatialDM will be pinned separately once the Python workflow is initialized.
