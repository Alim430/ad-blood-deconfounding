###############################################################################
# Step 22: pseudobulk (donor x cell-type) DE — the rigorous, reviewer-proof
# version of the within-cell-type DE (avoids single-cell pseudoreplication).
# Runs on Mac from the HVG export (2000 genes) + cells.csv (50 donors).
#   Rscript script/Paper1_Step22_pseudobulk.R
# Output: results/scRNA_pseudobulk_DE.csv
# NOTE: HVG export -> covers BTK/NFKBIA etc. but NOT IRAK1/RIPK3/LCK (not HVGs);
#       a full-gene version needs server aggregation from the 7GB object.
###############################################################################
suppressMessages({library(Matrix); library(DESeq2)})
setwd(".")

cat("Reading counts.mtx (genes x cells)...\n")
M <- as(readMM("results/ai/counts.mtx"), "CsparseMatrix")
genes <- readLines("results/ai/genes.txt"); rownames(M) <- genes
meta  <- read.csv("results/ai/cells.csv", stringsAsFactors=FALSE)
stopifnot(ncol(M)==nrow(meta))
cat(sprintf("  %d genes x %d cells\n", nrow(M), ncol(M)))

# donor = barcode prefix; pseudobulk group = donor x cell_type
meta$donor <- sub("_[ACGT]+-1$","",meta$barcode)
meta <- meta[!is.na(meta$condition) & meta$condition!="" , ]
M <- M[, meta$barcode %in% meta$barcode]                       # keep aligned (no-op safeguard)
meta$group <- paste(meta$donor, meta$cell_type, sep="__")
grp <- factor(meta$group)
ind <- sparse.model.matrix(~0+grp); colnames(ind) <- levels(grp) # cells x groups

cat("Aggregating to pseudobulk...\n")
pb <- as.matrix(M %*% ind); storage.mode(pb) <- "integer"        # genes x groups
gmeta <- unique(meta[,c("group","donor","cell_type","condition")])
rownames(gmeta) <- gmeta$group; gmeta <- gmeta[colnames(pb),]
gmeta$ncell <- as.integer(table(meta$group)[colnames(pb)])

# DE per cell type across donors (DESeq2), keeping pseudobulk samples with >=10 cells
res_all <- list()
for (ct in sort(unique(gmeta$cell_type))) {
  keep <- gmeta$cell_type==ct & gmeta$ncell>=10
  sm <- gmeta[keep,,drop=FALSE]
  sm$condition <- factor(sm$condition, levels=c("Control","AD"))
  if (length(unique(sm$condition))<2 || min(table(sm$condition))<3) {
    cat(sprintf("  [skip] %s: too few donors per group\n", ct)); next }
  cnt <- pb[, rownames(sm), drop=FALSE]; cnt <- cnt[rowSums(cnt)>=10,,drop=FALSE]
  dds <- DESeqDataSetFromMatrix(cnt, colData=sm, design=~condition)
  dds <- suppressMessages(DESeq(dds, quiet=TRUE))
  r <- as.data.frame(results(dds, contrast=c("condition","AD","Control")))
  r$gene <- rownames(r); r$cell_type <- ct
  r$n_AD <- sum(sm$condition=="AD"); r$n_Ctrl <- sum(sm$condition=="Control")
  res_all[[ct]] <- r
  cat(sprintf("  %-9s donors AD/Ctrl=%d/%d | padj<0.05: %d\n",
              ct, r$n_AD[1], r$n_Ctrl[1], sum(r$padj<0.05, na.rm=TRUE)))
}
out <- do.call(rbind, res_all)
write.csv(out, "results/scRNA_pseudobulk_DE.csv", row.names=FALSE)

# concordance with the cell-level FindMarkers DE
fm <- read.csv("results/scRNA_within_celltype_DE.csv")
sigPB <- subset(out, padj<0.05); sigFM <- fm   # fm already filtered to sig
key <- function(d) paste(d$cell_type, d$gene)
ov <- intersect(key(sigPB), key(sigFM))
cat(sprintf("\nPseudobulk sig (padj<0.05): %d gene-celltype pairs\n", nrow(sigPB)))
cat(sprintf("Confirmed vs cell-level DE (same cell type & gene): %d\n", length(ov)))
cat("\nKeystone genes at pseudobulk level:\n")
for (g in c("NFKBIA","BTK","ITGAM","TNFAIP3","PTGDS")) {
  sub <- out[out$gene==g, c("cell_type","log2FoldChange","padj")]
  if (nrow(sub)) { cat(" ",g,":\n"); print(sub[order(sub$padj),], row.names=FALSE) }
}
cat("\nSaved: results/scRNA_pseudobulk_DE.csv\n")
