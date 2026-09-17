# AGENTS.md

This repository contains the analysis code for a multi-omics fracture-healing study comparing Union and NonUnion outcomes.

The main goal is a transparent and reproducible scientific analysis. Prefer code that a researcher can read from top to bottom and understand without navigating through many layers of abstractions.

## General principles

1. Prefer simple, explicit analysis code over software-engineering abstractions.
2. Keep the biological and statistical logic visible in the analysis script.
3. Do not introduce infrastructure for hypothetical future use cases.
4. Follow the existing repository structure and naming conventions.
5. Keep each pull request focused on one concrete analysis step.
6. Preserve reproducibility without making the workflow unnecessarily complex.

## Repository structure

Executable analyses belong in:

- `scripts/massspec/`
- `scripts/olink/`
- `scripts/spatial/`
- `scripts/integration/`

Reusable functions may be placed in:

- `src/R/`
- `src/python/`

Only move code into `src/` when it is genuinely reused or when doing so clearly improves readability.

Do not create helper files merely to make an analysis script shorter.

Analysis numbering restarts within each modality, for example:

```text
scripts/spatial/
  01_import_qc.R
  02_...
```

Keep related analysis steps close together and easy to discover.

## Analysis code style

Prefer scripts that can be read sequentially:

```text
parameters
packages
input metadata
load data
annotate data
quality control
analysis
plots/tables
save outputs
```

Important scientific operations should remain visible in the main script. Examples include:

- filtering criteria
- sample inclusion/exclusion
- annotation matching
- construction of biological replicate IDs
- pseudobulk aggregation
- model formulas
- contrasts
- multiple-testing correction
- pathway-ranking statistics

Do not hide these operations behind chains of helper functions.

Small local functions are appropriate for repetitive technical tasks such as checking required columns, formatting labels, safe summary statistics, or repeated plotting code.

For one-off logic, explicit code is preferred.

Avoid deeply nested `if`/`else` blocks. Prefer simple validation followed by early `stop()` calls where appropriate.

Do not add defensive code for scenarios that do not occur in this dataset unless there is a clear reason.

## Configuration

Do not introduce YAML or other configuration files unless several scripts genuinely need to share the same parameters.

For a small number of analysis parameters, define them clearly near the top of the relevant script.

Sample-specific information belongs in metadata tables rather than config files.

## Paths and data

Never hard-code machine-specific absolute paths.

Use one of:

- project-relative paths
- explicit command-line arguments
- a single environment variable for an external data root

Do not commit raw sequencing data, Space Ranger outputs, large DIA-NN output files, generated Seurat objects, or other large intermediate files.

Generated outputs should go into the existing `results/`, `figures/`, or ignored processed-data directories.

## Reproducibility

R dependencies are managed with `renv`. Do not introduce another R dependency-management system.

When adding a dependency:

1. confirm that it is actually required;
2. add it to the existing environment;
3. avoid packages used only to replace simple base-R functionality.

Record software/session information for major analyses where useful.

Do not claim an analysis has been successfully executed if the required raw data were not available for testing.

## Statistical principles

The biological replicate, not an individual measurement unit, determines the level of statistical inference.

Do not use technical observations as independent biological replicates.

Examples:

- spatial transcriptomics: animal, not spot
- technical replicates: not independent animals
- repeated ROI/tissue measurements: remain linked to the originating animal

Keep exploratory and confirmatory analyses clearly distinguishable.

Report effect sizes and uncertainty where appropriate, not only p-values.

Multiple-testing correction should be used for genome/proteome-wide testing.

Do not silently change the statistical model, contrast, filtering rule, or biological unit from an existing analysis.

## Spatial transcriptomics

### Biological design

Standard Visium spots are not biological replicates.

For replicated inference:

- whole-section pseudobulk: one biological profile per animal
- ROI analysis: animal × ROI pseudobulk
- tissue analysis: animal × tissue pseudobulk

Use repeated-measures/blocking methods where multiple measurements from the same animal enter the same model.

### Exploratory time course

The D3/D7/D14/D21 screen with one Union and one NonUnion sample per time point is descriptive and hypothesis-generating.

Do not perform inferential Union-vs-NonUnion testing from spots in this exploratory dataset.

### Replicated D3/D14 analysis

D3 and D14 replicated analyses use biological replicates for inference.

Prioritize per-day Union vs NonUnion effects, effect sizes, pathway-level interpretation, and animal-level visualization.

Do not rely on spot-level pseudo-replication.

### ROI and tissue

ROI and tissue are separate spatial annotations.

Do not assume ROI defines tissue, tissue defines ROI, or that one is nested within the other.

Analyze them as separate spatial stratifications unless an explicit ROI × tissue analysis is justified by the data.

### Spatial annotations

Each sample/capture has a spot-level `spa.csv` in its sample folder.

Rows correspond to spatial barcodes/spots. Annotation matching must be performed by barcode within the corresponding sample/capture.

Keep sample-level metadata separate from spot-level annotations.

### Counts

Preserve raw counts.

Do not replace raw counts with SCT, integrated, Harmony-corrected, or other transformed values for pseudobulk differential-expression testing.

Normalization or transformed assays may be used later for visualization when appropriate.

### Communication analyses

CellChat, ligand-receptor screening, and SpatialDM answer related but different questions:

- CellChat: cell-type/pathway communication networks
- ligand-receptor screening: molecular ligand-receptor effects
- SpatialDM: spatially coordinated ligand-receptor expression

Run these at the animal level where possible and summarize across biological replicates.

Spatial ligand-receptor association is not proof of signaling.

## Mass spectrometry

Stay close to the standard `prolfqua` workflow where practical.

DIA-NN quantitative output is the primary mass-spec input.

Do not unnecessarily reimplement functionality already handled clearly by `prolfqua`.

Keep mass-spec preprocessing, QC, and differential analysis separate from Olink analysis.

## Olink

Analyze Olink independently from mass spectrometry before integration.

Because the Olink panel is small, do not perform standalone pathway enrichment as though it were a genome-wide assay.

Combined biological interpretation with mass-spec results can be performed later in `scripts/integration/`.

## Multi-omics integration

Do not mix modalities prematurely.

First produce interpretable results independently for spatial transcriptomics, mass spectrometry, and Olink.

Integration should operate on already-defined results or appropriately matched sample-level quantities.

Keep serum–spatial integration in `scripts/integration/`, not `scripts/spatial/`.

## Pull requests

Each PR should answer one concrete need.

Prefer a small number of changed files, one clearly described analysis step, readable diffs, explicit assumptions, and instructions for local execution.

Avoid PRs that simultaneously introduce infrastructure and several downstream analyses.

Before opening a PR:

1. inspect the existing code and follow its conventions;
2. check that paths are portable;
3. check that biological replicate IDs are correct;
4. check that raw counts are preserved where required;
5. run syntax/static checks;
6. run the analysis when the required data are available;
7. state clearly when full execution was not possible.

## Refactoring

Refactor only when it improves scientific readability or removes meaningful duplication.

Do not refactor working analysis code solely to make it more abstract, generic, or object-oriented.

When choosing between several helper functions that hide the workflow and a somewhat longer but readable sequential analysis script, prefer the readable sequential script.

Scientific auditability is more important than minimizing line count.
