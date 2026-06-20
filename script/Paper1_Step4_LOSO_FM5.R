###############################################################################
# Paper 1 — Step 4: LOSO Cross-Validation + FM5 Analytical Variability
#
# Section 1: LOSO (Leave-One-Study-Out) Cross-Validation
#             For each held-out study, re-run meta-analysis on remaining studies,
#             then check if the held-out study validates the pooled signature.
#             Reports per-gene stability and overall AUC.
#
# Section 2: FM5 Analytical Variability
#             Run DESeq2 vs edgeR vs limma-voom on GSE270454 (RNA-seq).
#             Genes with discordant results across methods get FM5 flag.
#
# Section 3: Update RS with LOSO stability dimension
#
# Section 4: Final consolidated biomarker table
#
# Prerequisites: results/03_meta_RS_FM.RData (from Step 3)
###############################################################################

setwd(".")

library(metafor)
library(limma)
library(data.table)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

load("results/03_meta_RS_FM.RData")
cat("Loaded Step 3 workspace.\n")
cat("Studies in effects_all:", paste(unique(effects_all$study), collapse = ", "), "\n")
cat("Genes in meta_results:", nrow(meta_results), "\n\n")


###############################################################################
# SECTION 1: LOSO Cross-Validation
#
# For each study S:
#   1. Remove S from effect size table
#   2. Re-run meta-analysis on remaining studies (genes in ≥2 remaining)
#   3. Select top N genes from leave-out meta
#   4. Check direction concordance in held-out study S
#   5. If possible, compute AUC using held-out expression data
#
# Output: per-gene LOSO stability score (fraction of leave-outs where
#         gene remains significant and direction-consistent)
###############################################################################

cat("=== SECTION 1: LOSO Cross-Validation ===\n\n")

all_studies <- unique(effects_all$study)
n_studies_total <- length(all_studies)
cat("Total studies:", n_studies_total, "\n")
cat("Studies:", paste(all_studies, collapse = ", "), "\n\n")

# Only meaningful if ≥3 studies (need ≥2 remaining after leave-out)
if (n_studies_total < 3) {
  cat("WARNING: Only", n_studies_total, "studies — LOSO needs ≥3.\n")
  cat("LOSO will use leave-one-out but remaining meta needs ≥2 studies per gene.\n")
  cat("Many genes will be dropped in each fold. This is expected.\n\n")
}

# Store per-gene LOSO results
loso_genes <- unique(effects_all$gene)
loso_matrix <- matrix(NA, nrow = length(loso_genes), ncol = n_studies_total)
rownames(loso_matrix) <- loso_genes
colnames(loso_matrix) <- all_studies

# For each leave-out fold:
#   loso_matrix[gene, study] = 1 if gene validates in held-out study
#                             = 0 if gene fails validation
#                             = NA if gene not testable in that fold

cat("Running LOSO folds...\n\n")

for (s_idx in seq_along(all_studies)) {
  held_out <- all_studies[s_idx]
  cat("--- Fold", s_idx, ": Hold out", held_out, "---\n")

  # Split data
  train_effects <- effects_all[effects_all$study != held_out, ]
  test_effects <- effects_all[effects_all$study == held_out, ]

  cat("  Training: ", nrow(train_effects), "gene-study pairs from",
      length(unique(train_effects$study)), "studies\n")
  cat("  Testing:  ", nrow(test_effects), "gene-study pairs from", held_out, "\n")

  # Count genes in ≥2 training studies
  train_gene_counts <- table(train_effects$gene)
  trainable_genes <- names(train_gene_counts[train_gene_counts >= 2])
  cat("  Genes with ≥2 training studies:", length(trainable_genes), "\n")

  if (length(trainable_genes) == 0) {
    cat("  No genes with ≥2 training studies — skipping fold.\n\n")
    next
  }

  # Re-run meta-analysis on training data
  train_meta <- data.frame(
    gene = character(),
    pooled_logFC = numeric(),
    pooled_pval = numeric(),
    stringsAsFactors = FALSE
  )

  for (g in trainable_genes) {
    dat <- train_effects[train_effects$gene == g, ]
    vi <- dat$SE^2

    tryCatch({
      fit <- rma(yi = dat$logFC, vi = vi, method = "REML")
      train_meta <- rbind(train_meta, data.frame(
        gene = g,
        pooled_logFC = fit$beta[1],
        pooled_pval = fit$pval,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      tryCatch({
        fit <- rma(yi = dat$logFC, vi = vi, method = "FE")
        train_meta <<- rbind(train_meta, data.frame(
          gene = g,
          pooled_logFC = fit$beta[1],
          pooled_pval = fit$pval,
          stringsAsFactors = FALSE
        ))
      }, error = function(e2) {})
    })
  }

  # FDR on training meta
  train_meta$adj_pval <- p.adjust(train_meta$pooled_pval, method = "BH")

  # Select significant genes from training
  sig_train <- train_meta[train_meta$adj_pval < 0.05, ]
  cat("  Significant in training meta:", nrow(sig_train), "\n")

  # Validate in held-out study
  test_genes <- intersect(sig_train$gene, test_effects$gene)
  cat("  Testable in held-out:", length(test_genes), "\n")

  if (length(test_genes) > 0) {
    for (g in test_genes) {
      train_direction <- sign(sig_train$pooled_logFC[sig_train$gene == g])
      test_row <- test_effects[test_effects$gene == g, ]

      # Validation criteria:
      # 1. Same direction as training meta
      # 2. Nominally significant (p < 0.05) in held-out study
      same_dir <- sign(test_row$logFC[1]) == train_direction
      nom_sig <- test_row$pvalue[1] < 0.05

      if (same_dir && nom_sig) {
        loso_matrix[g, held_out] <- 1  # validated
      } else if (same_dir) {
        loso_matrix[g, held_out] <- 0.5  # direction OK but not significant
      } else {
        loso_matrix[g, held_out] <- 0  # failed
      }
    }

    validated <- sum(loso_matrix[test_genes, held_out] == 1, na.rm = TRUE)
    partial <- sum(loso_matrix[test_genes, held_out] == 0.5, na.rm = TRUE)
    failed <- sum(loso_matrix[test_genes, held_out] == 0, na.rm = TRUE)
    cat("  Validated:", validated, "| Partial:", partial, "| Failed:", failed, "\n")
  }
  cat("\n")
}

# Compute per-gene LOSO stability score
# = mean score across all folds where gene was testable
loso_stability <- data.frame(
  gene = rownames(loso_matrix),
  loso_score = apply(loso_matrix, 1, function(x) mean(x, na.rm = TRUE)),
  loso_folds_tested = apply(loso_matrix, 1, function(x) sum(!is.na(x))),
  loso_folds_validated = apply(loso_matrix, 1, function(x) sum(x == 1, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

# Remove genes never tested
loso_stability <- loso_stability[loso_stability$loso_folds_tested > 0, ]
loso_stability <- loso_stability[order(-loso_stability$loso_score), ]

cat("\nLOSO Stability Summary:\n")
cat("Genes tested in ≥1 fold:", nrow(loso_stability), "\n")
cat("Genes validated in ALL folds:", sum(loso_stability$loso_score == 1, na.rm = TRUE), "\n")
cat("Mean LOSO score:", round(mean(loso_stability$loso_score, na.rm = TRUE), 3), "\n")

cat("\nTop 20 most stable genes:\n")
print(head(loso_stability, 20))

write.csv(loso_stability, "results/LOSO_stability.csv", row.names = FALSE)
cat("Saved: results/LOSO_stability.csv\n")


###############################################################################
# SECTION 2: FM5 — Analytical Variability
#
# Run three different DE methods on GSE270454 (RNA-seq):
#   1. DESeq2 (already done in Step 2)
#   2. edgeR (quasi-likelihood F-test)
#   3. limma-voom
#
# Compare direction and significance. Genes with discordant results
# across ≥2 methods → FM5 flag.
###############################################################################

cat("\n\n=== SECTION 2: FM5 Analytical Variability ===\n\n")

# Load RNA-seq counts
load("results/GSE270454_processed.RData")

# Set up comparison: AD vs NonDemented (ASM+ASO)
meta_270 <- meta_270454
meta_270$group <- ifelse(meta_270$condition == "AD", "AD",
                  ifelse(meta_270$condition %in% c("ASM", "ASO"), "NonDemented", "Other"))
keep_idx <- meta_270$group %in% c("AD", "NonDemented")
counts_fm5 <- counts_270454[, keep_idx]
meta_fm5 <- meta_270[keep_idx, ]
meta_fm5$group <- factor(meta_fm5$group, levels = c("NonDemented", "AD"))

# Remove zero-count genes
counts_fm5 <- counts_fm5[rowSums(counts_fm5) > 0, ]
cat("FM5 analysis: ", ncol(counts_fm5), "samples,", nrow(counts_fm5), "genes\n")
cat("Groups:", table(meta_fm5$group), "\n\n")

# --- Method 1: DESeq2 (already have res_270_ad_nd from Step 2) ---
cat("Method 1: DESeq2 (from Step 2)...\n")
deseq2_res <- res_270_ad_nd[, c("gene", "log2FoldChange", "pvalue", "padj")]
colnames(deseq2_res) <- c("gene", "logFC_deseq2", "pval_deseq2", "padj_deseq2")
cat("  Genes:", nrow(deseq2_res), "\n")

# --- Method 2: edgeR ---
cat("Method 2: edgeR...\n")
if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR", ask = FALSE)
}
library(edgeR)

dge <- DGEList(counts = counts_fm5, group = meta_fm5$group)
keep_edger <- filterByExpr(dge)
dge <- dge[keep_edger, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design_edger <- model.matrix(~ group, data = meta_fm5)
dge <- estimateDisp(dge, design_edger)

# Quasi-likelihood F-test (recommended for small samples)
fit_edger <- glmQLFit(dge, design_edger)
qlf <- glmQLFTest(fit_edger, coef = 2)  # AD vs NonDemented
edger_res <- topTags(qlf, n = Inf)$table
edger_res$gene <- rownames(edger_res)
edger_res <- edger_res[, c("gene", "logFC", "PValue", "FDR")]
colnames(edger_res) <- c("gene", "logFC_edger", "pval_edger", "padj_edger")
cat("  Genes:", nrow(edger_res), "\n")

# --- Method 3: limma-voom ---
cat("Method 3: limma-voom...\n")
dge_voom <- DGEList(counts = counts_fm5, group = meta_fm5$group)
keep_voom <- filterByExpr(dge_voom)
dge_voom <- dge_voom[keep_voom, , keep.lib.sizes = FALSE]
dge_voom <- calcNormFactors(dge_voom)

design_voom <- model.matrix(~ group, data = meta_fm5)
v <- voom(dge_voom, design_voom, plot = FALSE)
fit_voom <- lmFit(v, design_voom)
fit_voom <- eBayes(fit_voom)
voom_res <- topTable(fit_voom, coef = 2, number = Inf)
voom_res$gene <- rownames(voom_res)
voom_res <- voom_res[, c("gene", "logFC", "P.Value", "adj.P.Val")]
colnames(voom_res) <- c("gene", "logFC_voom", "pval_voom", "padj_voom")
cat("  Genes:", nrow(voom_res), "\n")

# --- Merge all three methods ---
cat("\nMerging results across methods...\n")
fm5_merged <- merge(deseq2_res, edger_res, by = "gene", all = TRUE)
fm5_merged <- merge(fm5_merged, voom_res, by = "gene", all = TRUE)
cat("Genes in all three methods:", sum(complete.cases(fm5_merged)), "\n")

# Direction concordance
fm5_merged$dir_deseq2 <- sign(fm5_merged$logFC_deseq2)
fm5_merged$dir_edger <- sign(fm5_merged$logFC_edger)
fm5_merged$dir_voom <- sign(fm5_merged$logFC_voom)

# Count how many methods agree on direction
fm5_merged$n_agree <- apply(fm5_merged[, c("dir_deseq2", "dir_edger", "dir_voom")], 1,
                             function(x) {
                               x <- x[!is.na(x)]
                               if (length(x) < 2) return(NA)
                               max(table(x))
                             })

fm5_merged$all_agree <- fm5_merged$n_agree == 3

# FM5 flag: direction discordance OR significance discordance
# Significance: significant in one method but not in others (at nominal p < 0.05)
fm5_merged$sig_deseq2 <- fm5_merged$pval_deseq2 < 0.05
fm5_merged$sig_edger <- fm5_merged$pval_edger < 0.05
fm5_merged$sig_voom <- fm5_merged$pval_voom < 0.05

fm5_merged$n_sig <- rowSums(fm5_merged[, c("sig_deseq2", "sig_edger", "sig_voom")],
                             na.rm = TRUE)

# FM5 = direction discordance (not all agree) OR mixed significance (1 of 3 sig)
fm5_merged$FM5 <- (!fm5_merged$all_agree) | (fm5_merged$n_sig == 1)
fm5_merged$FM5[is.na(fm5_merged$FM5)] <- FALSE

cat("\nFM5 Summary:\n")
cat("Direction concordance (all 3 agree):",
    sum(fm5_merged$all_agree, na.rm = TRUE), "of",
    sum(!is.na(fm5_merged$all_agree)), "\n")
cat("FM5 flagged:", sum(fm5_merged$FM5), "\n")

write.csv(fm5_merged, "results/FM5_analytical_variability.csv", row.names = FALSE)
cat("Saved: results/FM5_analytical_variability.csv\n")


###############################################################################
# SECTION 3: Update RS with LOSO + FM5
#
# Add D8 (LOSO stability) to RS computation.
# Update FM classification with FM5.
###############################################################################

cat("\n\n=== SECTION 3: Updated RS + FM ===\n\n")

# Merge LOSO scores into rs_data
rs_updated <- rs_data

# Add LOSO score (D8)
loso_match <- match(rs_updated$gene, loso_stability$gene)
rs_updated$D8_loso <- loso_stability$loso_score[loso_match]
rs_updated$D8_loso[is.na(rs_updated$D8_loso)] <- 0  # not tested = 0

# Add FM5
fm5_match <- match(rs_updated$gene, fm5_merged$gene)
rs_updated$FM5_analytical <- FALSE
rs_updated$FM5_analytical[!is.na(fm5_match)] <- fm5_merged$FM5[fm5_match[!is.na(fm5_match)]]

# Recompute RS with 8 dimensions
dim_cols_v2 <- c("D1_consistency", "D2_effect_size", "D3_stat_strength",
                  "D4_homogeneity", "D5_cross_platform", "D6_replication",
                  "D7_stage_spec", "D8_loso")

dim_matrix_v2 <- rs_updated[, dim_cols_v2]
for (col in dim_cols_v2) {
  na_idx <- is.na(dim_matrix_v2[[col]])
  if (any(na_idx)) {
    dim_matrix_v2[[col]][na_idx] <- median(dim_matrix_v2[[col]], na.rm = TRUE)
  }
  rng <- range(dim_matrix_v2[[col]], na.rm = TRUE)
  if (rng[2] > rng[1]) {
    dim_matrix_v2[[col]] <- (dim_matrix_v2[[col]] - rng[1]) / (rng[2] - rng[1])
  }
}

pca_v2 <- prcomp(dim_matrix_v2, center = TRUE, scale. = TRUE)
cat("Updated PCA variance explained:\n")
var_exp_v2 <- pca_v2$sdev^2 / sum(pca_v2$sdev^2) * 100
print(round(var_exp_v2[1:min(5, length(var_exp_v2))], 1))

cat("\nPC1 loadings:\n")
print(round(pca_v2$rotation[, 1], 3))

pc1_v2 <- pca_v2$x[, 1]
if (cor(pc1_v2, dim_matrix_v2$D3_stat_strength) < 0) pc1_v2 <- -pc1_v2
rs_updated$RS_v2 <- (pc1_v2 - min(pc1_v2)) / (max(pc1_v2) - min(pc1_v2)) * 100

# Update primary FM with FM5
for (i in seq_len(nrow(rs_updated))) {
  if (rs_updated$FM5_analytical[i] && rs_updated$primary_FM[i] == "Robust") {
    rs_updated$primary_FM[i] <- "FM5_Analytical"
  }
}

cat("\nUpdated FM Distribution:\n")
print(table(rs_updated$primary_FM))

cat("\nRS_v2 distribution:\n")
cat("  Mean:", round(mean(rs_updated$RS_v2), 1), "\n")
cat("  Median:", round(median(rs_updated$RS_v2), 1), "\n")

write.csv(rs_updated, "results/RS_FM_final.csv", row.names = FALSE)
cat("Saved: results/RS_FM_final.csv\n")


###############################################################################
# SECTION 4: Final Consolidated Biomarker Table
###############################################################################

cat("\n\n=== SECTION 4: Final Biomarker Table ===\n\n")

# Build the definitive table for the paper
final_table <- rs_updated[, c("gene", "n_studies", "pooled_logFC", "pooled_pval",
                                "adj_pval", "I2", "RS_v2",
                                "D1_consistency", "D2_effect_size", "D3_stat_strength",
                                "D4_homogeneity", "D5_cross_platform", "D6_replication",
                                "D7_stage_spec", "D8_loso",
                                "direction_consistent", "primary_FM",
                                "FM1_platform", "FM2_population", "FM3_stage",
                                "FM4_temporal", "FM5_analytical")]
colnames(final_table)[colnames(final_table) == "RS_v2"] <- "RS"
final_table <- final_table[order(-final_table$RS), ]

# Tier classification
final_table$tier <- "Tier3_Exploratory"
final_table$tier[final_table$RS >= 50 & final_table$adj_pval < 0.05] <- "Tier2_Candidate"
final_table$tier[final_table$RS >= 70 & final_table$adj_pval < 0.05 &
                   final_table$direction_consistent] <- "Tier1_Robust"

cat("Tier distribution:\n")
print(table(final_table$tier))

cat("\nTier 1 (Robust) biomarkers:\n")
t1 <- final_table[final_table$tier == "Tier1_Robust", ]
if (nrow(t1) > 0) {
  print(t1[, c("gene", "RS", "pooled_logFC", "adj_pval", "n_studies", "I2", "primary_FM")])
} else {
  cat("  None at RS ≥ 70. Showing top candidates (RS ≥ 60):\n")
  t1_relaxed <- final_table[final_table$RS >= 60 & final_table$adj_pval < 0.05 &
                              final_table$direction_consistent, ]
  if (nrow(t1_relaxed) > 0) {
    print(head(t1_relaxed[, c("gene", "RS", "pooled_logFC", "adj_pval",
                               "n_studies", "I2", "primary_FM")], 30))
  }
}

cat("\nTop 50 genes by RS (regardless of tier):\n")
print(head(final_table[, c("gene", "RS", "tier", "pooled_logFC", "adj_pval",
                            "n_studies", "I2", "primary_FM")], 50))

write.csv(final_table, "results/Paper1_final_biomarker_table.csv", row.names = FALSE)
cat("\nSaved: results/Paper1_final_biomarker_table.csv\n")

# Save workspace
save.image("results/04_LOSO_FM5_final.RData")
cat("Saved: results/04_LOSO_FM5_final.RData\n")
