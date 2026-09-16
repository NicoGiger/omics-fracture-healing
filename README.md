# omics-fracture-healing

Paper-facing analysis repository for the fracture-healing study, integrating spatial transcriptomics and serum proteomics.

## Repository structure

- `scripts/` — executable analysis entry points. Assay-specific workflows are grouped by modality and numbered independently from `01_...` within each subdirectory.
  - `scripts/massspec/` — DIA-NN / mass-spectrometry workflow.
  - `scripts/olink/` — Olink workflow.
  - `scripts/spatial/` — spatial-transcriptomics workflow.
  - `scripts/integration/` — cross-omics analyses.
  - `scripts/setup_R.R` and `scripts/prepare_serum_metadata.R` — shared setup/metadata entry points.
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

The spatial workflow is planned as a staged analysis that keeps exploratory and replicated evidence separate:

1. Exploratory D3/D7/D14/D21 time course (`n = 1` per condition/time point): whole-section overview, ROI-level description, and tissue-specific analysis.
2. Replicated D3/D14 whole-section pseudobulk analysis (`n = 3` per group/time point): sample-level PCA, differential expression, and pathway analysis.
3. Replicated ROI-level and tissue-specific pseudobulk/cell-composition analyses.
4. Cell-type deconvolution and spot-level pathway activity.
5. Spatial-organization analyses, including Moran's I.
6. Broad CellChat analysis of cell-cell communication.
7. Broad ligand-receptor screening.
8. Broad SpatialDM analysis for spatially coordinated ligand-receptor interactions.
9. Visium HD analysis for high-resolution localization and mechanistic exploration.

### Serum proteomics

Mass spectrometry and Olink are treated as separate assays through QC and differential analysis. They are integrated only downstream at the biological-interpretation stage.

#### Shared serum metadata

The supplied study metadata contains assay flags and the paired baseline/endpoint design. `scripts/prepare_serum_metadata.R` converts it to clean assay-specific tables with:

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
Rscript scripts/prepare_serum_metadata.R /path/to/Metadata.csv
```

This writes `metadata/massspec_samples.csv` and `metadata/olink_samples.csv`.

#### Mass spectrometry: DIA-NN -> prolfqua

The MS branch stays as close as possible to the standard [`prolfqua`](https://github.com/fgcz/prolfqua) workflow. Custom code is limited to DIA-NN input adaptation and study-specific design/contrast handling:

`DIA-NN -> AnalysisConfiguration -> setup_analysis() -> LFQData -> Transformer/Plotter -> Model -> Contrasts`

The stable `prolfqua` v1.5.0 release is pinned in `renv`.

The original DIA-NN report is large and should remain local. `scripts/massspec/01_extract_diann_quant.R` reads only the columns needed to compare available protein-level quantities (`PG.MaxLFQ`, `PG.Normalised`, `PG.Quantity`) and filtering information:

```bash
Rscript scripts/massspec/01_extract_diann_quant.R /path/to/report.tsv data/proteomics/massspec/diann_extract
```

It writes a compact `Run x Protein.Group` table plus field/consistency summaries. The current checked wide matrix was generated from DIA-NN `PG.Normalised`; the canonical protein quantity should be fixed only after inspecting the local DIA-NN extract.

MS analysis:

```bash
Rscript scripts/massspec/02_qc.R
# inspect results/proteomics/massspec/model_matrix_columns.txt
# add named contrasts to config/massspec.yml
Rscript scripts/massspec/03_differential.R
```

Defaults are conservative: log2 transformation of linear DIA-NN abundances, no extra normalization, no imputation, and no complete-case requirement. Process-QC samples are retained for QC but excluded from differential modelling. Repeated animals use a prolfqua random-intercept mixed model by default. The fixed effects use a cell-means parameterization (`0 + analysis_cell`) so day-specific paired changes and between-condition contrasts are explicit rather than hidden in a large interaction formula.

#### Olink Target 48 Mouse quantified concentrations

The supplied Olink workbook is a quantified-concentration export in `pg/mL`, not an NPX export. It therefore has a separate parser and statistical workflow rather than being forced through prolfqua.

Olink analysis:

```bash
Rscript scripts/olink/01_qc.R
# inspect results/proteomics/olink/sample_design_counts.csv
# inspect results/proteomics/olink/model_matrix_columns.txt
# add named contrasts to config/olink.yml
Rscript scripts/olink/02_differential.R
```

The scaffold parses the workbook's assay annotation rows (`Assay`, `Uniprot ID`, `OlinkID`, `Unit`), sample-level `QC Warning`, and quantified concentrations. Concentrations are log2 transformed for modelling. Missing/censored values are not imputed by default and no arbitrary assay-missingness threshold is imposed until the assay-specific missingness pattern is reviewed.

Differential analysis uses `limma` with the same explicit `analysis_cell` parameterization as MS. Repeated baseline/endpoint samples are handled with `duplicateCorrelation` and `animal_id` as the block. This also allows the observed Olink design to remain estimable when some condition x endpoint cells are absent.

#### Setup

Copy `config/paths.example.yml` to `config/paths.yml`, edit local paths, and initialize the R environment once:

```bash
Rscript scripts/setup_R.R
```

### Integration

Cross-omics analyses are performed only after assay-specific analyses, using `animal_id` as the shared experimental-unit key where measurements are matched. Serum-to-spatial analyses belong in `scripts/integration/` rather than either assay-specific directory.

## Reproducibility

R dependencies are managed with `renv`. Python dependencies for SpatialDM will be pinned separately once the Python workflow is initialized.
