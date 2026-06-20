# Data access

**No raw data is redistributed in this repository.** All datasets are public; ADNI is
controlled-access. Below are the accessions and how to obtain each.

## Bulk blood transcriptomes (discovery)
| Dataset | Accession | Platform | Samples | How to get |
|---|---|---|---|---|
| AddNeuroMed 1 | **GSE63060** | Illumina HumanHT-12 v3 | 329 (AD/MCI/CTL) | GEO (open) — `GEOquery::getGEO("GSE63060")` |
| AddNeuroMed 2 | **GSE63061** | Illumina HumanHT-12 v4 | 388 (AD/MCI/CTL) | GEO (open) |
| Third cohort  | **GSE140829** | Illumina (RNA) | 587 | GEO (open) |

## Independent validation
| Dataset | Access | Notes |
|---|---|---|
| **ADNI** gene expression (Affymetrix U219), n=343 | **Controlled-access** | Apply at https://adni.loni.usc.edu → Data & Samples → request access; sign the Data Use Agreement. **Do not redistribute.** This repo contains only derived summary tables (e.g., per-gene effect sizes), never ADNI-level data. |

## Single-cell RNA-seq + paired TCR
| Dataset | Accession | Notes |
|---|---|---|
| PBMC scRNA-seq + TCR | **GSE226602** | GEO (open); 50 donors, 270,828 cells. PBMC → **no neutrophils** (relevant to single-cell caveats in the manuscript). |

## Reproduction order
1. Download GEO datasets via the `script/Paper1_*Load*`/`download_robust.sh` helpers (or `GEOquery`).
2. Apply for ADNI separately; place its expression matrix where `Paper1_Step13_ADNI_validation.R` expects it (path noted at the top of that script).
3. Run `script/Paper1_Step2_*` … `Paper1_Step33b_*` in order.

## R packages (key)
limma, metafor, DESeq2, Seurat, edgeR, sva; BisqueRNA, MuSiC, SingleCellExperiment, Biobase
(for the multi-method deconvolution robustness check, Step32). See `sessionInfo.txt` for exact versions.

## Python packages
See `requirements.txt`.
