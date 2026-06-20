###############################################################################
# Paper 1 — Step 7: Single-Cell PBMC Analysis (GSE226602)  [LAYER 2]  v2 = RDS
#
# Decides, for each whole-blood SURVIVOR (Step 6) and each plasma-cfRNA hit
# (Step 8), whether it is:
#   - CELL-TYPE-SPECIFIC (one lineage)         -> composition-driven
#   - DIFFERENTIALLY EXPRESSED WITHIN a cell type (Control vs AD) -> cell-INTRINSIC
# This is the core decomposition of the reframed paper.
#
# DATA (already downloaded, ~8 GB):
#   scRNA/GSE226602_rna_raw_counts.rds        (gene x cell sparse counts OR Seurat)
#   scRNA/GSE226602_rna_lognorm_expression.rds
#   (optional) scRNA/GSE226602_series_matrix.txt.gz  (GSM -> diagnosis key)
#
# RUN ON THE CLOUD SERVER (more RAM). Needs Seurat:
#   CONDA_SOLVER=classic conda install -n scrna -c conda-forge r-seurat r-singler r-celldex
#   Rscript script/Paper1_Step7_scRNA.R
#
# >>> EDIT the CONFIG block once you know how diagnosis is encoded (barcodes vs metadata).
###############################################################################

## ---------------- CONFIG (EDIT THESE) -------------------------------------
SC_DIR     <- "scRNA"                                   # dir with the .rds files
COUNTS_RDS <- file.path(SC_DIR, "GSE226602_rna_raw_counts.rds")
LOGN_RDS   <- file.path(SC_DIR, "GSE226602_rna_lognorm_expression.rds")
# how to get Control/AD per cell. Options resolved in order:
#   1) if the loaded object is a Seurat with a diagnosis column -> set COND_COL
#   2) else derive from the cell barcode prefix -> set COND_FROM_BARCODE regex+map
SAMPLE_COND_CSV   <- "scRNA/sample_condition.csv"  # columns: sample,condition (Control/AD)
BARCODE_SUFFIX    <- "_[ACGT]+-1$"   # strip from a cell barcode -> sample prefix (e.g. G1010_y2)
## --------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(Matrix)
  ok_seurat <- requireNamespace("Seurat", quietly = TRUE)
  if (ok_seurat) library(Seurat)
  library(ggplot2)
}))
stopifnot(ok_seurat)  # install r-seurat into the scrna env first
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

## ---------------- genes of interest (from bulk + plasma) ------------------
if (file.exists("results/goi_genes.txt")) {
  goi <- trimws(unlist(strsplit(readLines("results/goi_genes.txt"), "[ ,\t]+"))); goi <- goi[goi!=""]
} else {
  gs <- character(0); gp <- character(0)
  if (file.exists("results/biomarker_survival_after_adjustment.csv")) {
    s <- read.csv("results/biomarker_survival_after_adjustment.csv", stringsAsFactors=FALSE)
    gs <- s$gene[s$survives %in% c(TRUE,"TRUE","True") & s$tier_orig=="Tier1_Robust"]
  }
  if (file.exists("results/cfRNA_translation_table.csv")) {
    c8 <- read.csv("results/cfRNA_translation_table.csv", stringsAsFactors=FALSE)
    gp <- c8$gene[c8$plasma_translatable %in% c(TRUE,"TRUE","True")]
  }
  goi <- unique(c(gs, gp))
}
cat("Genes of interest:", length(goi), "\n")

## ---------------- load + build Seurat ------------------------------------
cat("Loading", COUNTS_RDS, "...\n")
obj_raw <- readRDS(COUNTS_RDS)

if (inherits(obj_raw, "Seurat")) {
  pbmc <- obj_raw
} else {
  pbmc <- CreateSeuratObject(counts = as(obj_raw, "dgCMatrix"),
                             min.cells = 3, min.features = 200)
  rm(obj_raw); gc()
  pbmc <- NormalizeData(pbmc, verbose = FALSE)   # log-normalize from raw counts
}
cat("Cells:", ncol(pbmc), " Genes:", nrow(pbmc), "\n")

## ---------------- attach diagnosis (Control / AD) ------------------------
get_condition <- function(pbmc) {
  stopifnot(file.exists(SAMPLE_COND_CSV))
  sc <- read.csv(SAMPLE_COND_CSV, stringsAsFactors=FALSE)
  map <- setNames(sc$condition, sc$sample)
  pre <- sub(BARCODE_SUFFIX, "", colnames(pbmc))    # cell barcode -> sample prefix
  unname(map[pre])
}
pbmc$condition <- get_condition(pbmc)
cat("Condition counts:\n"); print(table(pbmc$condition, useNA="ifany"))
if (all(is.na(pbmc$condition)))
  stop("Diagnosis not resolved — inspect head(colnames(pbmc)) and set CONFIG (COND_COL or AD/CT patterns).")

## ---------------- QC + clustering ----------------------------------------
CKPT <- "results/_pbmc_clustered.rds"
if (file.exists(CKPT)) {
  cat("Resuming from clustering checkpoint (skipping QC/cluster/UMAP).\n")
  pbmc <- readRDS(CKPT)
} else {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern="^MT-")
  pbmc <- subset(pbmc, nFeature_RNA>200 & nFeature_RNA<5000 & percent.mt<15)
  # (already log-normalized at load)
  pbmc <- FindVariableFeatures(pbmc, nfeatures=2000, verbose=FALSE)
  pbmc <- ScaleData(pbmc, verbose=FALSE)
  pbmc <- RunPCA(pbmc, npcs=30, verbose=FALSE)
  pbmc <- FindNeighbors(pbmc, dims=1:20, verbose=FALSE)
  pbmc <- FindClusters(pbmc, resolution=0.5, verbose=FALSE)
  pbmc <- tryCatch(RunUMAP(pbmc, dims=1:20, verbose=FALSE),
                   error=function(e){cat("UMAP skipped (",conditionMessage(e),")\n"); pbmc})
  saveRDS(tryCatch(DietSeurat(pbmc, dimreducs=intersect(c("pca","umap"), Reductions(pbmc))),
                   error=function(e) pbmc), CKPT)
  cat("Saved clustering checkpoint -> re-runs skip the slow steps.\n")
}

## ---------------- cell-type annotation -----------------------------------
annotated <- FALSE
if (requireNamespace("SingleR",quietly=TRUE) && requireNamespace("celldex",quietly=TRUE)) {
  ref <- tryCatch(celldex::MonacoImmuneData(), error=function(e) NULL)
  if (!is.null(ref)) {
    pred <- SingleR::SingleR(test=GetAssayData(pbmc,"data"), ref=ref, labels=ref$label.main)
    pbmc$cell_type <- pred$labels; annotated <- TRUE; cat("Annotated via SingleR/Monaco.\n")
  }
}
if (!annotated) {
  mk <- list(
    CD4T=c("CD3D","CD3E","IL7R","CD4"), CD8T=c("CD8A","CD8B","GZMK"),
    NK=c("NKG7","GNLY","KLRD1","NCAM1"), Bcell=c("CD79A","MS4A1","CD19","BANK1"),
    Monocyte=c("CD14","LYZ","FCN1","S100A8"), DC=c("FCER1A","CLEC10A","CLEC9A"),
    Neutrophil=c("FCGR3B","CSF3R","CXCR2"), Platelet=c("PPBP","PF4"))
  avg <- AverageExpression(pbmc, group.by="seurat_clusters", assays="RNA")$RNA
  colnames(avg) <- sub("^g","",colnames(avg))   # AverageExpression prepends 'g' to numeric cluster ids
  call <- sapply(colnames(avg), function(cl){
    sco <- sapply(mk, function(g) mean(avg[intersect(g, rownames(avg)), cl], na.rm=TRUE))
    if (all(is.na(sco))) "Unknown" else names(which.max(sco)) })
  ct <- unname(call[as.character(pbmc$seurat_clusters)]); ct[is.na(ct)] <- "Unknown"
  pbmc$cell_type <- ct; cat("Annotated via markers.\n")
}
cat("Cell types:\n"); print(table(pbmc$cell_type))
saveRDS(pbmc, "results/GSE226602_seurat_processed.rds")
if ("umap" %in% names(pbmc@reductions)) try(ggsave("figures/Fig_scRNA_UMAP.png",
       DimPlot(pbmc, group.by="cell_type", label=TRUE, repel=TRUE)+ggtitle("GSE226602 PBMC"),
       width=7, height=6, dpi=150), silent=TRUE)

## ---------------- Q1: cell-type proportions Control vs AD -----------------
prop <- as.data.frame.matrix(table(pbmc$cell_type, pbmc$condition))
prop$cell_type <- rownames(prop)
write.csv(prop, "results/scRNA_celltype_counts.csv", row.names=FALSE)
cat("\nQ1 cell-type counts by condition:\n"); print(prop)

## ---------------- Q2: cell-type specificity of GOI -----------------------
goi_in <- intersect(goi, rownames(pbmc))
cat("\nGOI present in scRNA:", length(goi_in), "/", length(goi), "\n")
spec <- NULL
if (length(goi_in)>0) {
  avg <- AverageExpression(pbmc, features=goi_in, group.by="cell_type", assays="RNA")$RNA
  spec <- data.frame(gene=rownames(avg),
                     top_celltype=colnames(avg)[max.col(avg)],
                     specificity=apply(avg,1,function(x) max(x)/(sum(x)+1e-9)))
  write.csv(spec, "results/scRNA_goi_specificity.csv", row.names=FALSE)
  ggsave("figures/Fig_scRNA_goi_dotplot.png",
         DotPlot(pbmc, features=head(goi_in,30), group.by="cell_type")+RotatedAxis(),
         width=11, height=6, dpi=150)
}

## ---------------- Q3: within-cell-type Control-vs-AD DE ------------------
intrinsic <- list()
Idents(pbmc) <- "cell_type"
for (ct in unique(pbmc$cell_type)) {
  sub <- subset(pbmc, idents=ct)
  if (length(unique(na.omit(sub$condition)))<2) next
  if (sum(table(sub$condition)>=3)<2) next
  Idents(sub) <- "condition"
  m <- tryCatch(FindMarkers(sub, ident.1="AD", ident.2="Control",
            features=intersect(goi,rownames(sub)), logfc.threshold=0, min.pct=0.05, verbose=FALSE),
            error=function(e) NULL)
  if (!is.null(m) && nrow(m)>0){ m$cell_type<-ct; m$gene<-rownames(m); intrinsic[[ct]]<-m }
}
if (length(intrinsic)>0){
  idf <- do.call(rbind,intrinsic); idf <- idf[order(idf$p_val_adj),]
  write.csv(idf,"results/scRNA_within_celltype_DE.csv",row.names=FALSE)
  cat("\nQ3 cell-INTRINSIC DE (GOI still DE within a cell type):\n")
  print(head(idf[,c("gene","cell_type","avg_log2FC","p_val_adj")],25))
  cat(sprintf("\n%d GOI are significant in preliminary cell-level tests (padj<0.05 within a lineage).\n",
              length(unique(idf$gene[idf$p_val_adj<0.05]))))
}
save.image("results/07_scRNA.RData")
cat("\nSTEP 7 done. Interpretation: high specificity + NO within-CT DE => composition;",
    "Final inference requires donor-level pseudobulk; see Step22/Step22b/Step22c.\n")
