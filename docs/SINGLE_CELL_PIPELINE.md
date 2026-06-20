# Single-cell replication pipeline — AD blood scRNA pseudobulk meta-analysis

Goal: test whether the cell-intrinsic null (GSE226602: ~6 genes at donor-level) **replicates**
across independent AD blood scRNA cohorts, and whether the **143 bulk survivors** show any
concordant cell-intrinsic signal. Honest framing: GSE226602 (50 donors) is the best-powered
cohort; the others (5–6 donors) are direction-replication only, not power boosters.

## Stage 0 — datasets

| Accession | Donors (AD/C) | Cells | Type | Download |
|---|---|---|---|---|
| GSE226602 (primary, have it) | ~50 | 270k | PBMC+TCR | already processed |
| GSE181279 (Xu & Jia 2021) | 5 (3/2) | 36,849 | PBMC+TCR | GEO supplementary |
| Xiong 2021 (EMM) | 6 (4/2) | — | PBMC | check GEO/GSA accession in paper |
| GSE Gate 2020 (Nature) | small | — | PBMC | GEO supplementary |

Whole-blood (neutrophil-containing) AD scRNA: not available as of search — note as limitation.

## Stage 1 — download (prefer PROCESSED matrices; never FASTQ)

```bash
# example: GSE181279 supplementary (10x mtx or counts table)
mkdir -p scRNA_rep/GSE181279 && cd scRNA_rep/GSE181279
# list supplementary files:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181279  -> "Supplementary file"
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181279/suppl/GSE181279_RAW.tar"
tar -xvf GSE181279_RAW.tar
# read counts in R/Seurat (Read10X or read.table depending on format)
```
Robust-download tips: `wget -c` (resume), verify md5, prefer the `_RAW.tar` or a processed
`*_counts.csv.gz` over BAM/FASTQ. If only FASTQ exist, that dataset is not worth the cost here.

## Stage 2 — per-dataset preprocessing (server, Seurat) — reuse Step7 logic

For each cohort, produce a Seurat object with: `counts`, `cell_type`, `condition`, `donor`.
```r
library(Seurat)
cnt <- Read10X(...)                      # or read the counts matrix
obj <- CreateSeuratObject(cnt, min.cells=3, min.features=200)
obj[["pct.mt"]] <- PercentageFeatureSet(obj, "^MT-")
obj <- subset(obj, nFeature_RNA>200 & nFeature_RNA<6000 & pct.mt<15)
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA() |>
       FindNeighbors(dims=1:20) |> FindClusters(res=0.5)
# annotate by markers: CD14/LYZ=Monocyte, CD3D/CD3E=T(CD4/CD8 by CD4/CD8A),
#   MS4A1/CD79A=Bcell, NKG7/GNLY=NK, PPBP=Platelet
# donor + condition come from the GEO sample metadata (map barcode prefix -> donor -> Dx)
saveRDS(obj, "scRNA_rep/<GSE>_seurat.rds")
```

## Stage 3 — pseudobulk per dataset (donor × cell type) — reuse Step22b/c

```r
# aggregate raw counts to donor x cell_type, then DESeq2 AD vs Control per cell type
# (identical to script/ai/Step22b_server_pseudobulk.R + Step22c_mac_pseudobulk_DE.R)
# output: results/pseudobulk_<GSE>.csv  (gene, log2FoldChange, pvalue, padj, cell_type)
```

## Stage 4 — cross-dataset meta-analysis of the pseudobulk effects

For each cell type, meta-analyse the per-dataset pseudobulk log2FC with metafor (same engine
as the bulk meta), giving a pooled cell-intrinsic effect per gene per cell type.
```r
library(metafor)
# effects = rbind of per-dataset (gene, log2FC, lfcSE, cell_type, study)
# per gene+cell_type: rma(yi=log2FC, vi=lfcSE^2, method="REML")
# -> results/scRNA_pseudobulk_META.csv ; count padj<0.05 per cell type
```

## Stage 5 — the two questions to answer

1. **Does the null replicate?** Across cohorts, is the meta pseudobulk still ~null
   (few genes, dominated by Hb/Ig contamination)? If yes → the cell-intrinsic null is
   replicated and credible (not a one-cohort fluke).
2. **Do the 143 bulk survivors show cell-intrinsic signal?** Intersect the 143 with the meta
   pseudobulk; test enrichment of low p-values vs background. Expected: no — confirming the
   bulk survivors act via composition, not cell-intrinsic expression.

## Honest expectations & caveats
- GSE226602 (50 donors) dominates power; GSE181279/Xiong (5–6 donors) are direction checks only.
- A concordant near-null across cohorts is the likely (and publishable) outcome: "the peripheral
  AD cell-intrinsic signal is minimal and reproducibly so."
- This does NOT need GPU (pseudobulk + metafor are CPU); only Seurat preprocessing needs RAM.
- Neutrophils remain untestable (PBMC) — state as limitation; whole-blood scRNA is the real gap.
```
