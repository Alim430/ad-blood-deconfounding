###############################################################################
# Paper 1 — Step 2: Differential Expression Analysis
#
# Works with what we have NOW:
#   GSE140829: 587 samples (microarray, limma)
#   GSE270454: 45 samples (RNA-seq, DESeq2)
#
# GSE63060/63061 can be added later (Section 5) after browser download.
#
# Run in RStudio section by section.
###############################################################################

setwd(".")

raw_dir <- "./data/raw"

library(limma)
library(data.table)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Load processed data from Step 1
load("results/GSE140829_processed.RData")
load("results/GSE270454_processed.RData")

cat("GSE140829:", ncol(expr_140829_matched), "samples,",
    nrow(expr_140829_matched), "probes\n")
cat("GSE270454:", ncol(counts_270454), "samples,",
    nrow(counts_270454), "genes\n")


###############################################################################
# SECTION 1: GSE140829 — Preprocessing + limma DEG
###############################################################################

cat("\n=== SECTION 1: GSE140829 limma Analysis ===\n\n")

# 1a. Log2 transform (if not already log scale)
cat("Checking if data is log scale...\n")
cat("Range of values:", range(expr_140829_matched[1:100, 1:5], na.rm = TRUE), "\n")
max_val <- max(expr_140829_matched[1:1000, 1:10], na.rm = TRUE)
if (max_val > 100) {
  cat("Values > 100 detected — applying log2 transform\n")
  expr_log2 <- log2(expr_140829_matched + 1)
} else {
  cat("Data appears already log-transformed\n")
  expr_log2 <- expr_140829_matched
}
cat("Log2 range:", round(range(expr_log2[1:100, 1:5], na.rm = TRUE), 2), "\n")

# 1b. Filter low-expression probes
# Keep probes detected in at least 10% of samples
cat("\nFiltering low-expression probes...\n")
median_expr <- apply(expr_log2, 1, median, na.rm = TRUE)
keep <- median_expr > quantile(median_expr, 0.25)  # Drop bottom 25%
expr_filtered <- expr_log2[keep, ]
cat("Probes after filtering:", nrow(expr_filtered), "(removed",
    sum(!keep), "low-expression probes)\n")

# 1c. Quantile normalization
cat("Applying quantile normalization...\n")
expr_norm <- normalizeBetweenArrays(expr_filtered, method = "quantile")

# 1d. Collapse probes to genes (take max-expressing probe per gene)
cat("\nCollapsing probes to genes...\n")
probe_in_data <- probe_gene_map[probe_gene_map$probe_id %in% rownames(expr_norm), ]
cat("Probes with gene symbols in filtered data:", nrow(probe_in_data), "\n")

# For each gene, keep the probe with highest median expression
probe_medians <- data.frame(
  probe_id = probe_in_data$probe_id,
  symbol = probe_in_data$symbol,
  median_expr = apply(expr_norm[probe_in_data$probe_id, ], 1, median, na.rm = TRUE)
)
best_probes <- probe_medians[order(probe_medians$symbol, -probe_medians$median_expr), ]
best_probes <- best_probes[!duplicated(best_probes$symbol), ]
cat("Unique genes:", nrow(best_probes), "\n")

expr_genes <- expr_norm[best_probes$probe_id, ]
rownames(expr_genes) <- best_probes$symbol

# 1e. Set up design matrix — AD vs Control
cat("\nSetting up limma design matrix...\n")
meta <- meta_140829_matched

# Primary comparison: AD vs Control
# Also keep MCI for secondary analysis
meta$group <- factor(meta$diagnosis, levels = c("Control", "MCI", "AD"))
cat("Group distribution:\n")
print(table(meta$group))

# Design matrix
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

# Add covariates if available
if ("Sex" %in% colnames(meta) && "age_at_draw" %in% colnames(meta)) {
  meta$age_num <- as.numeric(meta$age_at_draw)
  meta$sex_bin <- ifelse(meta$Sex == "Male", 1, 0)
  design <- model.matrix(~ 0 + group + age_num + sex_bin, data = meta)
  colnames(design)[1:3] <- levels(meta$group)
  cat("Design includes: group + age + sex\n")
}

# 1f. Fit limma model
cat("\nFitting limma model...\n")
fit <- lmFit(expr_genes, design)

# Contrasts: AD vs Control, MCI vs Control, AD vs MCI
contrast_matrix <- makeContrasts(
  AD_vs_Control = AD - Control,
  MCI_vs_Control = MCI - Control,
  AD_vs_MCI = AD - MCI,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 1g. Extract results
cat("\n--- AD vs Control ---\n")
deg_AD_ctrl <- topTable(fit2, coef = "AD_vs_Control", number = Inf, sort.by = "P")
deg_AD_ctrl$gene <- rownames(deg_AD_ctrl)
cat("Total genes tested:", nrow(deg_AD_ctrl), "\n")
cat("Significant (adj.P < 0.05):", sum(deg_AD_ctrl$adj.P.Val < 0.05), "\n")
cat("  Upregulated:", sum(deg_AD_ctrl$adj.P.Val < 0.05 & deg_AD_ctrl$logFC > 0), "\n")
cat("  Downregulated:", sum(deg_AD_ctrl$adj.P.Val < 0.05 & deg_AD_ctrl$logFC < 0), "\n")
cat("\nTop 20 DEGs:\n")
print(head(deg_AD_ctrl[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 20))

cat("\n--- MCI vs Control ---\n")
deg_MCI_ctrl <- topTable(fit2, coef = "MCI_vs_Control", number = Inf, sort.by = "P")
deg_MCI_ctrl$gene <- rownames(deg_MCI_ctrl)
cat("Significant (adj.P < 0.05):", sum(deg_MCI_ctrl$adj.P.Val < 0.05), "\n")

cat("\n--- AD vs MCI ---\n")
deg_AD_MCI <- topTable(fit2, coef = "AD_vs_MCI", number = Inf, sort.by = "P")
deg_AD_MCI$gene <- rownames(deg_AD_MCI)
cat("Significant (adj.P < 0.05):", sum(deg_AD_MCI$adj.P.Val < 0.05), "\n")

# Save DEG results
write.csv(deg_AD_ctrl, "results/DEG_GSE140829_AD_vs_Control.csv", row.names = FALSE)
write.csv(deg_MCI_ctrl, "results/DEG_GSE140829_MCI_vs_Control.csv", row.names = FALSE)
write.csv(deg_AD_MCI, "results/DEG_GSE140829_AD_vs_MCI.csv", row.names = FALSE)
cat("\nSaved DEG results to results/\n")


###############################################################################
# SECTION 2: GSE270454 — DESeq2 Analysis
###############################################################################

cat("\n\n=== SECTION 2: GSE270454 DESeq2 Analysis ===\n\n")

# Install DESeq2 if not available
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  cat("Installing DESeq2...\n")
  BiocManager::install("DESeq2", ask = FALSE)
}
library(DESeq2)

# For this dataset: AD vs non-AD (ASM + ASO as controls, since no healthy controls)
# Or: AD vs MCI as a secondary comparison
meta_270 <- meta_270454

# Option 1: AD vs MCI (cleanest comparison)
cat("Running AD vs MCI comparison...\n")
ad_mci_idx <- meta_270$condition %in% c("AD", "MCI")
counts_sub <- counts_270454[, ad_mci_idx]
meta_sub <- meta_270[ad_mci_idx, ]
meta_sub$condition <- factor(meta_sub$condition, levels = c("MCI", "AD"))

# Remove genes with 0 counts
nonzero <- rowSums(counts_sub) > 0
counts_sub <- counts_sub[nonzero, ]
cat("Genes with nonzero counts:", nrow(counts_sub), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = meta_sub,
  design = ~ condition
)

dds <- DESeq(dds)
res_270_ad_mci <- results(dds, contrast = c("condition", "AD", "MCI"))
res_270_ad_mci <- as.data.frame(res_270_ad_mci)
res_270_ad_mci$gene <- rownames(res_270_ad_mci)
res_270_ad_mci <- res_270_ad_mci[order(res_270_ad_mci$pvalue), ]

cat("\nAD vs MCI results:\n")
cat("Total genes:", nrow(res_270_ad_mci), "\n")
cat("Significant (padj < 0.05):", sum(res_270_ad_mci$padj < 0.05, na.rm = TRUE), "\n")
cat("  Up:", sum(res_270_ad_mci$padj < 0.05 & res_270_ad_mci$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Down:", sum(res_270_ad_mci$padj < 0.05 & res_270_ad_mci$log2FoldChange < 0, na.rm = TRUE), "\n")

# Option 2: AD vs ASM+ASO (treating them as non-demented)
cat("\nRunning AD vs non-demented (ASM+ASO) comparison...\n")
meta_270$group2 <- ifelse(meta_270$condition == "AD", "AD",
                   ifelse(meta_270$condition %in% c("ASM", "ASO"), "NonDemented", "MCI"))
ad_nd_idx <- meta_270$group2 %in% c("AD", "NonDemented")
counts_sub2 <- counts_270454[, ad_nd_idx]
meta_sub2 <- meta_270[ad_nd_idx, ]
meta_sub2$group2 <- factor(meta_sub2$group2, levels = c("NonDemented", "AD"))

counts_sub2 <- counts_sub2[rowSums(counts_sub2) > 0, ]

dds2 <- DESeqDataSetFromMatrix(
  countData = counts_sub2,
  colData = meta_sub2,
  design = ~ group2
)
dds2 <- DESeq(dds2)
res_270_ad_nd <- as.data.frame(results(dds2, contrast = c("group2", "AD", "NonDemented")))
res_270_ad_nd$gene <- rownames(res_270_ad_nd)
res_270_ad_nd <- res_270_ad_nd[order(res_270_ad_nd$pvalue), ]

cat("AD vs NonDemented results:\n")
cat("Total genes:", nrow(res_270_ad_nd), "\n")
cat("Significant (padj < 0.05):", sum(res_270_ad_nd$padj < 0.05, na.rm = TRUE), "\n")

# Save
write.csv(res_270_ad_mci, "results/DEG_GSE270454_AD_vs_MCI.csv", row.names = FALSE)
write.csv(res_270_ad_nd, "results/DEG_GSE270454_AD_vs_NonDemented.csv", row.names = FALSE)
cat("Saved DEG results.\n")


###############################################################################
# SECTION 3: Cross-Platform Concordance Check
#
# Compare GSE140829 (microarray) vs GSE270454 (RNA-seq) DEG directions
# This is your cross-platform consistency metric for RS scoring
###############################################################################

cat("\n\n=== SECTION 3: Cross-Platform Concordance ===\n\n")

# Get genes present in both datasets
genes_140829 <- deg_AD_ctrl$gene  # Already has gene symbols
genes_270454 <- res_270_ad_nd$gene  # Gene symbols from count matrix

common_genes <- intersect(genes_140829, genes_270454)
cat("Genes in GSE140829:", length(genes_140829), "\n")
cat("Genes in GSE270454:", length(genes_270454), "\n")
cat("Common genes:", length(common_genes), "\n")

# Merge by gene
cross <- merge(
  deg_AD_ctrl[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  res_270_ad_nd[, c("gene", "log2FoldChange", "pvalue", "padj")],
  by = "gene", suffixes = c("_micro", "_rnaseq")
)
cat("Merged rows:", nrow(cross), "\n")

# Direction concordance
cross$same_direction <- sign(cross$logFC) == sign(cross$log2FoldChange)
cat("\nDirection concordance (all genes):",
    round(mean(cross$same_direction, na.rm = TRUE) * 100, 1), "%\n")

# Concordance among significant genes
sig_both <- cross$adj.P.Val < 0.05 & cross$padj < 0.05
if (sum(sig_both, na.rm = TRUE) > 0) {
  cat("Direction concordance (significant in both):",
      round(mean(cross$same_direction[sig_both], na.rm = TRUE) * 100, 1), "%\n")
  cat("Genes significant in both:", sum(sig_both, na.rm = TRUE), "\n")
}

# Concordance among top 500 in each
top500_micro <- head(deg_AD_ctrl$gene[order(deg_AD_ctrl$P.Value)], 500)
top500_rnaseq <- head(res_270_ad_nd$gene[order(res_270_ad_nd$pvalue)], 500)
top_overlap <- intersect(top500_micro, top500_rnaseq)
cat("\nTop 500 overlap:", length(top_overlap), "genes\n")

if (length(top_overlap) > 0) {
  top_cross <- cross[cross$gene %in% top_overlap, ]
  cat("Direction concordance (top 500 overlap):",
      round(mean(top_cross$same_direction, na.rm = TRUE) * 100, 1), "%\n")
}

# Correlation plot data
write.csv(cross, "results/cross_platform_concordance.csv", row.names = FALSE)
cat("\nSaved: results/cross_platform_concordance.csv\n")


###############################################################################
# SECTION 4: Volcano Plots + Cross-Platform Scatter
###############################################################################

cat("\n\n=== SECTION 4: Figures ===\n\n")

# 4a. Volcano plot — GSE140829 AD vs Control
cat("Generating volcano plot...\n")
png("figures/Fig1_Volcano_GSE140829_AD_vs_Control.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 2))

plot(deg_AD_ctrl$logFC, -log10(deg_AD_ctrl$P.Value),
     pch = 20, cex = 0.5, col = "grey60",
     xlab = "log2 Fold Change (AD vs Control)",
     ylab = "-log10(P-value)",
     main = "GSE140829: AD vs Control (587 samples)",
     xlim = c(-2, 2))

# Highlight significant genes
sig <- deg_AD_ctrl$adj.P.Val < 0.05 & abs(deg_AD_ctrl$logFC) > 0.3
points(deg_AD_ctrl$logFC[sig], -log10(deg_AD_ctrl$P.Value[sig]),
       pch = 20, cex = 0.7,
       col = ifelse(deg_AD_ctrl$logFC[sig] > 0, "red", "blue"))

# Label top genes
top_genes <- head(deg_AD_ctrl[sig, ], 15)
if (nrow(top_genes) > 0) {
  text(top_genes$logFC, -log10(top_genes$P.Value),
       labels = top_genes$gene, cex = 0.7, pos = 4, col = "black")
}

abline(h = -log10(0.05), lty = 2, col = "grey40")
abline(v = c(-0.3, 0.3), lty = 2, col = "grey40")
legend("topright",
       legend = c(paste("Up:", sum(sig & deg_AD_ctrl$logFC > 0)),
                  paste("Down:", sum(sig & deg_AD_ctrl$logFC < 0)),
                  "NS"),
       col = c("red", "blue", "grey60"), pch = 20, cex = 0.9)
dev.off()
cat("Saved: figures/Fig1_Volcano_GSE140829_AD_vs_Control.png\n")

# 4b. Cross-platform scatter
if (nrow(cross) > 0) {
  cat("Generating cross-platform scatter...\n")
  png("figures/Fig2_CrossPlatform_Scatter.png", width = 700, height = 700)
  par(mar = c(5, 5, 4, 2))

  plot(cross$logFC, cross$log2FoldChange,
       pch = 20, cex = 0.4, col = "grey50",
       xlab = "logFC GSE140829 (Microarray)",
       ylab = "log2FC GSE270454 (RNA-seq)",
       main = paste("Cross-Platform Concordance:", length(common_genes), "genes"))

  # Highlight concordant significant
  if (any(sig_both, na.rm = TRUE)) {
    points(cross$logFC[sig_both], cross$log2FoldChange[sig_both],
           pch = 20, cex = 0.8, col = "red")
  }

  abline(h = 0, v = 0, lty = 2, col = "grey40")
  abline(lm(log2FoldChange ~ logFC, data = cross), col = "blue", lwd = 2)

  r <- cor(cross$logFC, cross$log2FoldChange, use = "complete.obs")
  legend("topleft", legend = paste("r =", round(r, 3)), bty = "n", cex = 1.2)
  dev.off()
  cat("Saved: figures/Fig2_CrossPlatform_Scatter.png\n")
}


###############################################################################
# SECTION 5: Load GSE63060/GSE63061 (after Safari download)
#
# Run this AFTER you download the files in Safari to data/raw/
###############################################################################

cat("\n\n=== SECTION 5: GSE63060/GSE63061 (run after browser download) ===\n\n")

# Check if files exist and are complete
for (gse in c("GSE63060", "GSE63061")) {
  f <- file.path(raw_dir, paste0(gse, "_series_matrix.txt.gz"))
  if (file.exists(f)) {
    fsize <- file.info(f)$size
    if (fsize > 50e6) {
      cat(gse, ": READY (", round(fsize/1e6, 1), "MB)\n")

      tryCatch({
        gse_data <- getGEO(filename = f)
        expr_data <- exprs(gse_data)
        meta_data <- pData(gse_data)
        cat("  Expression:", nrow(expr_data), "probes x", ncol(expr_data), "samples\n")

        ch1_cols <- grep(":ch1$", colnames(meta_data), value = TRUE)
        for (col in ch1_cols) {
          cat("  ", col, ":", paste(head(unique(meta_data[[col]]), 5), collapse = " | "), "\n")
        }

        # Save processed data
        save(expr_data, meta_data,
             file = paste0("results/", gse, "_processed.RData"))
        cat("  Saved: results/", gse, "_processed.RData\n")

      }, error = function(e) {
        cat("  Error loading:", e$message, "\n")
      })
    } else {
      cat(gse, ": File too small (", round(fsize/1e6, 1), "MB) — still truncated\n")
    }
  } else {
    cat(gse, ": Not downloaded yet\n")
    cat("  Download:", paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", gse, "nnn/", gse, "/matrix/", gse, "_series_matrix.txt.gz\n"))
  }
}


###############################################################################
# SECTION 6: Save Full Workspace
###############################################################################

cat("\n\n=== Saving Everything ===\n")
save.image("results/02_DEG_analysis.RData")
cat("Saved: results/02_DEG_analysis.RData\n")

cat("\n================================================================\n")
cat("  ANALYSIS COMPLETE — NEXT STEPS\n")
cat("================================================================\n\n")
cat("1. Check the volcano plot: figures/Fig1_Volcano_GSE140829_AD_vs_Control.png\n")
cat("2. Check cross-platform scatter: figures/Fig2_CrossPlatform_Scatter.png\n")
cat("3. Review DEG lists in results/ folder\n")
cat("4. Download GSE63060/63061 in Safari, then re-run Section 5\n")
cat("5. Next script: Meta-analysis + RS scoring (Paper1_Step3_MetaAnalysis.R)\n")
