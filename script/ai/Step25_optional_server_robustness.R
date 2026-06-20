###############################################################################
# Step 25 [SERVER, OPTIONAL]: remaining single-cell robustness checks.
# None expected to change the conclusion; for reviewer completeness only.
# Needs the 90GB instance (full Seurat object). CPU is enough (no GPU).
#   Rscript script/ai/Step25_optional_server_robustness.R
###############################################################################
setwd(Sys.getenv("AD_PAPER_ROOT", "."))
suppressWarnings(suppressMessages({library(Seurat); library(Matrix)}))
obj <- readRDS("results/GSE226602_seurat_processed.rds")
md <- obj@meta.data
cat("metadata columns:", paste(colnames(md), collapse=", "), "\n\n")

## 1. Batch / run metadata present? and association with condition --------------
bcol <- grep("batch|run|lane|sample|orig.ident|pool", colnames(md), value=TRUE, ignore.case=TRUE)
cat("candidate batch columns:", paste(bcol, collapse=", "), "\n")
for (b in bcol) {
  t <- try(table(md[[b]], md$condition), silent=TRUE)
  if (!inherits(t,"try-error") && nrow(t)>1 && nrow(t)<60) {
    p <- tryCatch(chisq.test(t)$p.value, error=function(e) NA)
    cat(sprintf("  %s vs condition: chisq p=%.2g (low p = batch-condition association)\n", b, p)) }
}

## 2. Doublet rate (scDblFinder if available) ----------------------------------
if (requireNamespace("scDblFinder", quietly=TRUE)) {
  suppressMessages({library(scDblFinder); library(SingleCellExperiment)})
  sce <- as.SingleCellExperiment(DietSeurat(obj, assays="RNA"))
  sce <- scDblFinder(sce)
  cat(sprintf("\nDoublet rate: %.1f%% (%d/%d)\n",
      100*mean(sce$scDblFinder.class=="doublet"), sum(sce$scDblFinder.class=="doublet"), ncol(sce)))
  cat("  by condition:\n"); print(round(100*prop.table(table(sce$scDblFinder.class, obj$condition),2),1))
} else cat("\n[skip] scDblFinder not installed (BiocManager::install('scDblFinder') to run)\n")

## 3. PCA-based donor-level AUC (light cross-check of scVI 0.52/0.68) -----------
suppressMessages({library(pROC)})
obj <- NormalizeData(obj, verbose=FALSE)
donor <- sub("_[ACGTN]+-1$","", colnames(obj))
# pseudobulk by donor (mean lognorm of HVGs)
obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)
X <- GetAssayData(obj, layer="data")[VariableFeatures(obj), ]
pb <- t(apply(X, 1, function(v) tapply(v, donor, mean)))   # genes x donor
cond <- tapply(as.character(obj$condition), donor, function(z) z[1])[colnames(pb)]
y <- as.integer(cond=="AD")
pc <- prcomp(t(pb), scale.=TRUE)$x[,1:10]
library(glmnet)
set.seed(1); fold <- sample(rep(1:5, length.out=length(y)))
pr <- numeric(length(y))
for (k in 1:5) { tr<-fold!=k
  cv <- cv.glmnet(pc[tr,], y[tr], family="binomial", alpha=0.5)
  pr[!tr] <- predict(cv, pc[!tr,], s="lambda.min", type="response") }
cat(sprintf("\nPCA donor-level AUC (10 PCs, 5-fold) = %.3f  (cf. scVI cell-level 0.52, donor pseudobulk 0.68)\n",
            as.numeric(pROC::auc(y, pr, quiet=TRUE))))
cat("\n[done] Step25 optional robustness. None of these should change the conclusion.\n")
