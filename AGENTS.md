# AGENTS.md

This repository contains analysis code for a multi-omics fracture-healing study comparing Union and NonUnion outcomes.

## Core principles

- Prefer simple, explicit analysis code over abstraction.
- Keep biological and statistical logic visible in the analysis script.
- Do not introduce infrastructure for hypothetical future use cases.
- Keep each pull request focused on one concrete analysis step.
- Prefer readable sequential scripts over chains of helper functions.
- Use `src/R/` or `src/python/` only for code that is genuinely reused or clearly improves readability.

## Analysis scripts

Executable analyses belong in `scripts/<modality>/` and numbering restarts within each modality.

A good analysis script should read roughly top-to-bottom as:

```text
parameters
packages
metadata/input
load data
annotate data
quality control
analysis
plots/tables
save outputs
```

Keep important scientific decisions visible in the main script, including filtering, sample inclusion, replicate IDs, pseudobulk construction, model formulas, contrasts, and multiple-testing correction.

Use small local helpers only for repetitive technical plumbing. Avoid deeply nested `if`/`else`; prefer simple checks and early `stop()` calls.

Do not add defensive code for scenarios that do not occur in this dataset without a clear reason.

## Configuration, paths, and dependencies

- Do not add YAML/config files unless several scripts genuinely share the same parameters.
- Put small numbers of analysis parameters near the top of the relevant script.
- Keep sample-specific information in metadata tables.
- Never hard-code machine-specific absolute paths.
- Use project-relative paths, a command-line argument, or one external data-root environment variable.
- Do not commit raw data, Space Ranger outputs, DIA-NN outputs, generated Seurat objects, or other large intermediates.
- R dependencies are managed with `renv`; do not introduce another dependency-management system.
- Do not add a package when simple base R is sufficient.
- Do not claim an analysis was run successfully if the required data were unavailable.

## Statistical principles

- The animal is the biological replicate unless explicitly documented otherwise.
- Never treat spots or technical replicates as independent biological replicates.
- Keep exploratory and confirmatory analyses distinct.
- Use multiple-testing correction for genome/proteome-wide testing.
- Do not silently change filtering rules, statistical models, contrasts, or the biological unit of inference.

## Spatial transcriptomics

- Preserve raw counts for pseudobulk differential-expression analyses.
- Whole-section inference: one pseudobulk profile per animal.
- ROI inference: animal × ROI pseudobulk.
- Tissue inference: animal × tissue pseudobulk.
- ROI and tissue are separate annotations; do not assume one is nested within the other.
- Each sample/capture has a spot-level `spa.csv`; match annotations to spots by barcode within that capture.
- Keep sample-level metadata separate from spot-level annotations.
- The D3/D7/D14/D21 Union/NonUnion time-course screen with n=1 per group/time is descriptive only.
- Replicated D3/D14 analyses use biological replicates for inference.
- Do not use SCT, Harmony/integration values, or other transformed values in place of raw counts for pseudobulk DE.
- Run spatial communication analyses per animal where possible and summarize across biological replicates.

## Proteomics and integration

- Keep mass spectrometry and Olink analyses separate before integration.
- Stay close to the standard `prolfqua` workflow for DIA-NN mass-spec analysis where practical.
- Do not treat the small Olink panel as a genome-wide assay for standalone enrichment.
- Put cross-modality analyses in `scripts/integration/`, not inside a modality-specific folder.

## Pull requests and refactoring

- Each PR should solve one concrete need with a small, readable diff.
- Follow existing repository conventions before adding new structure.
- Refactor only when it improves scientific readability or removes meaningful duplication.
- Do not create helper files merely to shorten an analysis script.
- If choosing between several helper layers and a somewhat longer but clear sequential script, prefer the clear script.

Scientific auditability is more important than minimizing line count.
