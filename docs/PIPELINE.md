# Reproducible Pipeline

This repository contains the public, cleaned pipeline for the v7 manuscript:

> Peripheral blood molecular signatures of Alzheimer's disease are predominantly explained by myeloid leukocyte composition.

The public package is deliberately conservative: raw data, ADNI files, generated result tables, large intermediate objects, detailed work logs, and historical exploratory modules are not redistributed. Scripts read data from local `data/` and write regenerated outputs to local `results/` and `figures/pub/`.

## Data

See `docs/data_access.md` for accession numbers and access notes.

- GEO bulk blood cohorts: GSE140829, GSE63060, GSE63061, GSE270454.
- GEO PBMC scRNA/TCR cohort: GSE226602.
- ADNI blood expression: controlled access through LONI; scripts require local ADNI files but no ADNI data are redistributed.

## Main Analysis Order

| Order | Script | Purpose |
|---|---|---|
| 1 | `Paper1_Load_Via_API.R`, `Paper1_Process_Available.R` | Load and harmonise available GEO data. |
| 2 | `Paper1_Step2_DEG_Analysis.R` | Initial per-cohort differential expression. |
| 3 | `Paper1_Step3_MetaAnalysis.R` | Random-effects meta-analysis nominating the Tier-1 candidate biomarkers. |
| 4 | `Paper1_Step4_LOSO_FM5.R` | Leave-one-study-out cross-validation used to characterise the candidate set. |
| 5 | `Paper1_Step6_Adjusted_Reanalysis.R` | Primary de-confounding analysis: batch/composition/age/sex adjustment and survivor attrition. |
| 6 | `Paper1_Step6b_LOCO_honest.R`, `Paper1_Step6c_uniform_deconfound.R` | Leak-free LOCO and uniform-model sensitivity checks. |
| 7 | `Paper1_Step9_CellComposition.R` | Reproducible myeloid-up/lymphoid-down composition shift. |
| 8 | `Paper1_Step13_ADNI_validation.R`, `Paper1_Step13b_ADNI_deconfound.R` | Independent ADNI replication with and without composition adjustment. |
| 9 | `Paper1_Step7_scRNA.R`, `Paper1_Step22_pseudobulk.R`, `script/ai/Step22b_server_pseudobulk.R`, `script/ai/Step22c_mac_pseudobulk_DE.R` | Single-cell processing and donor-level pseudobulk tests. Cell-level tests are retained only as the pseudoreplication contrast. |
| 10 | `script/ai/Step20a_export_h5.R`, `script/ai/Step20b_embedding.py` | scVI latent-space null check. |
| 11 | `Paper1_Step24_TCR.py` | Donor-level, depth-corrected TCR repertoire analysis. |
| 12 | `script/ai/Step17a_export_panel.R`, `script/ai/Step17b_ml_shap.py` | Interpretable ML ceiling and SHAP analysis in ADNI. |
| 13 | `Paper1_Step26_publication_figures.py`, `Paper1_Step27_remaining_figures.py`, `Paper1_Step29_multipanel.py`, `Paper1_Step30_graphical_abstract.py`, `Paper1_Step31_concept_paradox.py`, `Paper1_Step32_deconv_robustness.R`, `Paper1_Step33b_deconv_fig_real.py` | Publication figures and deconvolution robustness figure. |

## What Is Not Published

- Raw expression matrices and controlled-access ADNI files.
- Generated `results/` tables and large `.RData`/`.rds`/`.h5ad`/`.mtx` intermediates.
- The detailed internal work log.
- Historical full-pipeline scripts and exploratory controllability/virtual-KO/K562 perturbation modules that are not part of the verified v7 manuscript.

## Reproduction Note

Several later scripts intentionally depend on objects produced by earlier steps. If a file under `results/` is missing, run the preceding step rather than looking for it in the public repository.
