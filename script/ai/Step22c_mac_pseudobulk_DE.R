###############################################################################
# Step 22c [MAC]: full-gene pseudobulk DE from the aggregated object (DESeq2).
# Reads results/pseudobulk_fullgene.rds (downloaded from server).
#   Rscript script/ai/Step22c_mac_pseudobulk_DE.R
# Output: results/scRNA_pseudobulk_DE_fullgene.csv
###############################################################################
suppressWarnings(suppressMessages({library(DESeq2)}))
setwd(".")
x <- readRDS("results/pseudobulk_fullgene.rds"); pb <- x$pb; gm <- x$gm
cat(sprintf("%d genes x %d pseudobulk samples; cell types: %s\n",
            nrow(pb), ncol(pb), paste(sort(unique(gm$cell_type)),collapse=", ")))

res <- list()
for (ct in sort(unique(gm$cell_type))){
  keep <- gm$cell_type==ct & gm$ncell>=10; sm <- gm[keep,,drop=FALSE]
  sm$condition <- factor(sm$condition, levels=c("Control","AD"))
  if (length(unique(sm$condition))<2 || min(table(sm$condition))<3){ cat("  [skip]",ct,"\n"); next }
  cnt <- pb[,rownames(sm),drop=FALSE]; cnt <- cnt[rowSums(cnt)>=10,,drop=FALSE]
  dds <- DESeqDataSetFromMatrix(cnt, sm, ~condition)
  dds <- suppressMessages(DESeq(dds, quiet=TRUE))
  r <- as.data.frame(results(dds, contrast=c("condition","AD","Control")))
  r$gene<-rownames(r); r$cell_type<-ct; r$n_AD<-sum(sm$condition=="AD"); r$n_Ctrl<-sum(sm$condition=="Control")
  res[[ct]] <- r
  cat(sprintf("  %-9s AD/Ctrl=%d/%d  padj<0.05: %d\n", ct, r$n_AD[1], r$n_Ctrl[1], sum(r$padj<0.05,na.rm=TRUE)))
}
out <- do.call(rbind, res); write.csv(out, "results/scRNA_pseudobulk_DE_fullgene.csv", row.names=FALSE)
cat("\nKeystone genes (best padj across cell types):\n")
for (g in c("BTK","IRAK1","RIPK3","LCK","NFKBIA","ITGAM","IKBKG","PILRA","PTGDS")){
  s <- out[out$gene==g & !is.na(out$padj),]; s <- s[order(s$padj),]
  if (nrow(s)) cat(sprintf("  %-7s %-9s log2FC=%+.2f padj=%.3f\n", g, s$cell_type[1], s$log2FoldChange[1], s$padj[1]))
  else cat(sprintf("  %-7s not testable\n", g))
}
cat("\nSaved: results/scRNA_pseudobulk_DE_fullgene.csv\n")
