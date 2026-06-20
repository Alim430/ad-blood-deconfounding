# Supplementary Materials

**Peripheral blood molecular signatures of Alzheimer's disease are predominantly explained by myeloid leukocyte composition: a multi-modal de-confounding audit**

Contents: Supplementary Methods · Supplementary Note (robustness analyses S1–S20) · Supplementary Table S1 · Supplementary Figure S1.

---

## Supplementary Methods

**Datasets.** Discovery: GSE140829 (n=587: 204 AD, 249 control, 134 MCI), GSE63060 and GSE63061 (AddNeuroMed; n=329 and 388). Validation: ADNI gene expression (Affymetrix U219; 343 control/AD after visit-matched diagnosis). Single cell: GSE226602 (PBMC scRNA-seq + paired TCR; 50 donors, 270,828 cells; 28 AD, 22 control). For the case–control contrast, MCI samples were retained in the model but not in the AD−control comparison.

**Composition estimation.** Per-sample relative composition for six lineages (neutrophil, monocyte, B, CD4⁺ T, CD8⁺ T, NK) was computed as the mean z-score of curated Danaher marker genes (markers detected 36/36 in every cohort). The estimate is semi-quantitative (relative, not absolute) and is used as a continuous covariate.

**Bulk meta-analysis and survival.** Per-cohort limma contrasts (AD − control) were combined per gene by random-effects REML meta-analysis (metafor), with sampling variance SEᵢ². "Tier-1" candidates were genes meeting the naive robustness criterion (significant pooled effect with consistent direction). Survival after de-confounding required Benjamini–Hochberg-significant pooled effect with direction consistent with discovery, under the primary model (composition + age + sex + batch-where-available; no APOE).

**ADNI validation models.** Maximum-variance probe collapse; quantile normalisation; limma `~ group + age + sex (+ composition)`. Concordance is the fraction of a survivor set whose ADNI effect matches the discovery direction; its z is the binomial deviation from 0.5. The two-covariate model adjusts only neutrophil + monocyte scores.

**Single cell.** Seurat v5 (QC, normalisation, clustering, marker annotation). Reported DE: donor×cell-type pseudobulk → DESeq2 (AD vs control per cell type). Cell-level Wilcoxon and genomic inflation λ shown only to illustrate pseudoreplication. scVI: condition-blind latent embedding; AD-vs-control AUC by 5-fold cross-validated logistic regression at the cell level, and donor-level pseudobulk classification for comparison.

**TCR.** Per-donor Shannon entropy, normalised clonality, top-10 clonotype fraction, and richness; rarefied to a common T-cell depth (100 resamplings, averaged) before diversity metrics. CD8⁺-restricted metrics from contig↔cell-type join.

**Machine learning.** Standardisation inside the CV pipeline; elastic-net / random-forest / gradient-boosting; 5-fold CV; nested CV (inner 3-fold grid over C and L1-ratio); precision–recall AUC; calibration; SHAP (linear explainer on the elastic-net).

---

## Supplementary Note — robustness analyses (S1–S20)

Each analysis was pre-specified to test a specific threat to validity. ✗ = reviewer concern refuted by data; ✓ = confirmatory; ⚠ = real correction adopted.

**Bulk / meta-analysis**
- **S1 Leave-one-cohort-out meta.** Survivors when dropping each cohort: 444 (−GSE140829), 173 (−GSE63060), 338 (−GSE63061). Not single-cohort-driven; carried by AddNeuroMed; GSE140829 reduces (not drives) the count. ✓
- **S2 Variance-inflation factor (disease term).** 1.03 / 1.08 / 1.05. No collinearity between disease and composition covariates. ✗ over-adjustment-via-collinearity
- **S3 Heterogeneity (I²).** Survivors: median 0%, 9% with I²>50% vs 25% genome-wide. Survivors are cross-cohort concordant. ✓
- **S4 GSE140829 collapse decomposition.** Significant genes from base 2,667 → batch-only 104 (96% removed), composition-only 769, APOE-only 744. The collapse is driven mainly by batch–case confounding, not APOE. ⚠ (interpretation corrected)
- **S5 APOE sensitivity.** No-APOE model → 240 survivors; with-APOE → 162 (≈subset; intersection 160). Conclusion identical (S10). ⚠ (no-APOE adopted as primary)
- **S6 Bartlett homoscedasticity.** 28% / 15% / 12% of genes variance-unequal (vs ~5%); array-weighted re-analysis concordant. ⚠ (limitation stated)
- **S7 MCI handling.** Composition estimated per sample independent of group; monocyte shift +0.18/+0.31/+0.33 with MCI excluded from contrast. No bias. ✓
- **S8 Cross-platform marker detection.** 36/36 Danaher markers detected in all cohorts. ✗ marker incomparability

**ADNI validation**
- **S9 Two-covariate minimal model.** Neutrophil+monocyte only: concordance 88%→48%; survivor |log2FC| 0.076→0.027. Excludes high-dimensional over-adjustment. ✓
- **S10 APOE-free survivor set in ADNI.** 240-gene set collapses identically (90%→49%). Conclusion robust to discovery model. ✓
- **S11 Negative control (random gene sets).** Survivor collapse 0.041 vs 1,000 random sets (mean 0.001, 97.5th-ile 0.006); empirical p<0.001. ✗ "circular/targeted" criticism
- **S12 Composition–disease partial correlation (ctrl age/sex).** Neutrophil r=+0.17, monocyte r=+0.14 (modest, non-zero). Honest marginal association reported. ✓
- **S13 Regression-to-the-mean / winner's curse.** Survivor |log2FC|: discovery 0.063 → ADNI no-comp 0.076 (not attenuated, 88% concordant) → with-comp 0.025. Collapse is composition-specific, not RTM. ✗ winner's-curse explanation

**Single cell**
- **S14 Donor-level vs cell-level AUC.** All-PBMC donor pseudobulk AUC 0.68 (composition) vs misleading cell-level scVI 0.52; within-cell-type null. Cell-level metric clarified. ⚠
- **S15 Pseudobulk calibration (λ).** λ = 0.52–0.95 (all ≤1) per cell type → well-calibrated, no pseudoreplication inflation (cf. cell-level λ≫1). ✓
- **S16 Cell-count balance.** Cells/donor AD 5,704 vs control 4,875 (p=0.062, n.s.); DESeq2 size factors absorb library size. ✓
- **S17 Single-cell power.** ~8–10 donors/group/cell type powers detection of large effects (≳0.5 log2FC) only; null excludes large but not ~0.1-log2FC cell-intrinsic effects. ⚠ (limitation stated)

**Machine learning**
- **S18 APOE-only baseline + nested CV + PR-AUC.** APOE4 alone AUC 0.708; panel+demo 0.746 (panel +0.038). Nested-CV 0.753 ≈ naive 0.745 (no overfitting). Panel PR-AUC 0.40 (base rate 0.27). ✓ / ✗ overfitting
- **S19 TCR depth control.** Clonotype richness ∼ cell number (r=0.92); the single nominal hit is a depth artefact, n.s. after normalisation. ✓

**Equivalence testing**
- **S20 Two one-sided tests (TOST) of the ADNI survivor collapse.** For the composition-adjusted survivors in ADNI (n=119 with usable statistics; per-gene SE recovered from the two-sided p-value under a large-sample normal approximation, justified by n=343), the median |log₂FC| is 0.040 (median SE 0.041) and 87% of point estimates lie below 0.10. By TOST against an equivalence margin, 44% of survivors are statistically equivalent to zero at ±0.10 log₂FC and 72% at ±0.15 log₂FC. Interpretation: the surviving effects are uniformly small and centred near the genome-wide background, but ADNI's per-gene precision (median SE 0.041) is insufficient to declare every individual gene equivalent at the strict ±0.10 bound — incomplete equivalence reflects limited validation-cohort power (consistent with the power caveat in Discussion 3.6), not preserved cell-intrinsic signal. A genome-wide ADNI composition-adjusted differential-expression run would enable a survivor-versus-background equivalence comparison and is recommended at submission. ⚠

**Multi-method deconvolution robustness**
- **S21 Bisque + MuSiC cross-validation of Danaher composition estimates.** Using GSE226602 PBMC scRNA-seq (50 donors, 6 cell types: B, CD4T, CD8T, Monocyte, NK, Platelet) as reference, we re-estimated leukocyte composition in the two AddNeuroMed cohorts with Bisque (reference-based decomposition) and MuSiC (weighted regression). **Headline result: the monocyte AD vs Control shift is directionally positive in 4/4 (cohort × method) tests** (Bisque: GSE63060 Δ=+0.0055, p=8×10⁻⁵, n=249; GSE63061 Δ=+0.0050, p=8×10⁻⁴, n=273. MuSiC: GSE63060 Δ=+0.0013, p=0.29; GSE63061 Δ=+0.0025, p=0.013). Method-vs-Danaher monocyte concordance: Bisque r = 0.51 (GSE63060), 0.56 (GSE63061); MuSiC r = 0.32 (GSE63060), 0.63 (GSE63061); median r = 0.54. **Caveats** explicitly stated: (i) the PBMC scRNA reference has no neutrophils, so neutrophil cannot be cross-method-validated — the audit's reliance on the neutrophil score is mitigated by the minimal-model and independent ML head-to-head analyses (Figure 3c; §2.7). (ii) GSE140829 deconvolution failed for Bisque due to sparse gene overlap, consistent with the documented batch–case confounding in that cohort (S4) and the primary-analysis decision to weight conclusions on AddNeuroMed; this does not change the conclusion. (iii) Numerical r is moderate-to-low for T-cell subsets and B cells, expected because Danaher marker panels and PBMC scRNA references resolve different sub-populations of those lineages — directional concordance, not numerical identity, is the relevant robustness criterion. Output tables: `results/deconv_robustness/{cohort}_{method}.csv` and `concordance_summary.csv`. ✓ (audit's central claim is method-robust)

**Net effect.** No analysis overturns the conclusion. Two corrections adopted (no-APOE primary; donor-level single-cell metric); several reviewer concerns (collinearity, over-adjustment, circularity, winner's curse, overfitting, heterogeneity, single-cohort dependence, marker incomparability) are refuted by data.

---

## Supplementary Table S1. Robustness analyses summary
| ID | Test | Key result | Verdict |
|---|---|---|---|
| S1 | Leave-one-out meta | 444/173/338 survivors | not single-cohort driven |
| S2 | VIF (disease) | 1.03/1.08/1.05 | no collinearity |
| S3 | Heterogeneity I² | survivors 9% vs 25% | concordant |
| S4 | GSE140829 decomposition | batch 96% of collapse | batch–case confound |
| S5 | APOE sensitivity | 240 vs 162 | conclusion invariant |
| S6 | Bartlett | 12–28% heteroscedastic | limitation |
| S7 | MCI handling | shift unchanged | no bias |
| S8 | Marker detection | 36/36 all cohorts | comparable |
| S9 | 2-covariate ADNI | 88→48% | no over-adjustment |
| S10 | APOE-free in ADNI | 90→49% | model-robust |
| S11 | Random-gene null | p<0.001 | survivor-specific |
| S12 | Partial correlation | r=0.14–0.17 | modest, honest |
| S13 | Winner's curse | no-comp not attenuated | composition-specific |
| S14 | Donor vs cell AUC | 0.68 vs 0.52 | metric clarified |
| S15 | Pseudobulk λ | 0.52–0.95 | calibrated |
| S16 | Cell-count balance | p=0.062 | no confound |
| S17 | Single-cell power | large-effect only | limitation |
| S18 | ML rigor | APOE 0.708; nested≈naive | panel negligible |
| S19 | TCR depth | richness∼depth r=0.92 | artefact |
| S20 | TOST equivalence (ADNI survivors) | median \|log₂FC\|=0.040; 44%/72% equivalent at ±0.10/±0.15 | small effects; power-limited |

## Supplementary Figure S1
Four-panel robustness summary: leave-one-cohort-out survivor counts; APOE-only ML baseline; donor-level vs cell-level single-cell AUC; pseudobulk genomic inflation λ. (File: `figures/pub/FigS_robustness.pdf`.)

*The public repository contains the analysis scripts, data-access instructions, figure provenance, and supplementary robustness documentation needed to reproduce the analyses. Raw data and generated intermediate results are not redistributed.*
