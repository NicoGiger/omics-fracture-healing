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

Mass spectrometry and Olink are treated as separate assays through QC and differential analysis. They are integrated only downstream at the biological-interpretation stage.

#### Shared serum metadata

The supplied study metadata contains assay flags and the paired baseline/endpoint design. `scripts/00_prepare_serum_metadata.R` converts it to clean assay-specific tables with:

- `sample_id`
- `animal_id`
- `condition`
- measured `day`
- each animal's `endpoint_day` cohort
- `visit` (`baseline` / `endpoint`)
- `analysis_cell` (for explicit cell-means contrasts)
- `sample_type` (`biological` / `process_qc`)

Example:

```bash
Rscript scripts/00_prepare_serum_metadata.R /path/to/Metadata.csv
```

This writes `metadata/massspec_samples.csv` and `metadata/olink_samples.csv`.

#### Mass spectrometry: DIA-NN -> prolfqua

The MS branch stays as close as possible to the standard [`prolfqua`](https://github.com/fgcz/prolfqua) workflow. Custom code is limited to DIA-NN input adaptation and study-specific design/contrast handling:

`DIA-NN -> AnalysisConfiguration -> setup_analysis() -> LFQData -> Transformer/Plotter -> Model -> Contrasts`

The stable `prolfqua` v1.5.0 release is pinned in `renv`.

The original DIA-NN report is large and should remain local. `scripts/00_extract_diann_quant.R` reads only the columns needed to compare available protein-level quantities (`PG.MaxLFQ`, `PG.Normalised`, `PG.Quantity`) and filtering information:

```bash
Rscript scripts/00_extract_diann_quant.R /path/to/report.tsv data/proteomics/massspec/diann_extract
```

It writes a compact `Run x Protein.Group` table plus field/consistency summaries. The current checked wide matrix was generated from DIA-NN `PG.Normalised`; the canonical protein quantity should be fixed only after inspecting the local DIA-NN extract.

MS analysis:

```bash
Rscript scripts/01_ms_qc.R
# inspect results/proteomics/massspec/model_matrix_columns.txt
# add named contrasts to config/massspec.yml
Rscript scripts/02_ms_differential.R
```

Defaults are conservative: log2 transformation of linear DIA-NN abundances, no extra normalization, no imputation, and no complete-case requirement. Process-QC samples are retained for QC but excluded from differential modelling. Repeated animals use a prolfqua random-intercept mixed model by default. The fixed effects use a cell-means parameterization (`0 + analysis_cell`) so day-specific paired changes and between-condition contrasts are explicit rather than hidden in a large interaction formula.

#### Olink Target 48 Mouse quantified concentrations

The supplied Olink workbook is a quantified-concentration export in `pg/mL`, not an NPX export. It therefore has a separate parser and statistical workflow rather than being forced through prolfqua.

Olink analysis:

```bash
Rscript scripts/03_olink_qc.R
# inspect results/proteomics/olink/sample_design_counts.csv
# inspect results/proteomics/olink/model_matrix_columns.txt
# add named contrasts to config/olink.yml
Rscript scripts/04_olink_differential.R
```

The scaffold parses the workbook's assay annotation rows (`Assay`, `Uniprot ID`, `OlinkID`, `Unit`), sample-level `QC Warning`, and quantified concentrations. Concentrations are log2 transformed for modelling. Missing/censored values are not imputed by default and no arbitrary assay-missingness threshold is imposed until the assay-specific missingness pattern is reviewed.

Differential analysis uses `limma` with the same explicit `analysis_cell` parameterization as MS. Repeated baseline/endpoint samples are handled with `duplicateCorrelation` and `animal_id` as the block. This also allows the observed Olink design to remain estimable when some condition x endpoint cells are absent.

#### Setup

Copy `config/paths.example.yml` to `config/paths.yml`, edit local paths, and initialize the R environment once:

```bash
Rscript scripts/00_setup_R.R
```

### Integration

Cross-omics analyses are performed only after assay-specific analyses, using `animal_id` as the shared experimental-unit key where measurements are matched.

## Reproducibility

R dependencies are managed with `renv`. Python dependencies for SpatialDM will be pinned separately once the Python workflow is initialized.
