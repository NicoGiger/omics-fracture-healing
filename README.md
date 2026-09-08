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

1. Assay-specific QC and preprocessing
2. Protein-level temporal and outcome analyses
3. Pathway-level interpretation

### Integration

Cross-omics analyses are performed only after assay-specific analyses, using `animal_id` as the shared experimental-unit key where measurements are matched.

## Reproducibility

R dependencies will be managed with `renv`. Python dependencies for SpatialDM will be pinned separately once the Python workflow is initialized.
