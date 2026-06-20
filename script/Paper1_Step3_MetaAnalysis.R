###############################################################################
# Paper 1 — Step 3: Meta-Analysis + RS Scoring + FM Classification
#
# This script:
#   Section 1: Fix & load GSE63060/63061 (getGPL=FALSE)
#   Section 2: Run limma DEG on GSE63060/63061
#   Section 3: Standardize effect sizes across all cohorts
#   Section 4: Per-gene random effects meta-analysis (metafor REML)
#   Section 5: RS (Robustness Scoring) — PCA-driven composite
#   Section 6: FM (Failure Mode Classification) — FM1-FM5
#   Section 7: Summary tables + figures
#
# Prerequisites:
#   - results/02_DEG_analysis.RData (from Step 2)
#   - data/raw/GSE63060_series_matrix.txt.gz (62.8 MB, downloaded via Safari)
#   - data/raw/GSE63061_series_matrix.txt.gz (63.7 MB, downloaded via Safari)
#
# Run in RStudio section by section.
###############################################################################

setwd(".")

library(limma)
library(data.table)
library(GEOquery)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

raw_dir <- "./data/raw"

# Load Step 2 workspace
load("results/02_DEG_analysis.RData")
cat("Loaded Step 2 workspace.\n")
cat("Available DEG results: GSE140829, GSE270454\n\n")


###############################################################################
# SECTION 1: Load GSE63060 & GSE63061
#
# KEY FIX: Use getGPL=FALSE to prevent the HTTP request that was failing.
# The files are already downloaded (62.8 MB and 63.7 MB = complete).
###############################################################################

cat("=== SECTION 1: Load GSE63060 & GSE63061 ===\n\n")

# --- GSE63060 (AddNeuroMed cohort 1) ---
f63060 <- file.path(raw_dir, "GSE63060_series_matrix.txt.gz")
if (file.exists(f63060) && file.info(f63060)$size > 50e6) {
  cat("Loading GSE63060 (", round(file.info(f63060)$size/1e6, 1), "MB)...\n")

  gse63060 <- getGEO(filename = f63060, getGPL = FALSE)
  expr_63060 <- exprs(gse63060)
  meta_63060 <- pData(gse63060)

  cat("  Expression:", nrow(expr_63060), "probes x", ncol(expr_63060), "samples\n")

  # Show metadata fields
  ch1_cols_60 <- grep(":ch1$", colnames(meta_63060), value = TRUE)
  cat("  Metadata fields:\n")
  for (col in ch1_cols_60) {
    cat("    ", col, ":", paste(head(unique(meta_63060[[col]]), 5), collapse = " | "), "\n")
  }

  # Identify diagnosis column
  diag_col_60 <- NULL
  for (col in ch1_cols_60) {
    vals <- unique(tolower(meta_63060[[col]]))
    if (any(grepl("ad|alzheimer|control|mci|dementia", vals))) {
      diag_col_60 <- col
      cat("\n  Diagnosis column:", col, "\n")
      print(table(meta_63060[[col]]))
      break
    }
  }

  save(expr_63060, meta_63060, file = "results/GSE63060_processed.RData")
  cat("  Saved: results/GSE63060_processed.RData\n")
} else {
  cat("GSE63060 not available — skipping.\n")
}

# --- GSE63061 (AddNeuroMed cohort 2) ---
f63061 <- file.path(raw_dir, "GSE63061_series_matrix.txt.gz")
if (file.exists(f63061) && file.info(f63061)$size > 50e6) {
  cat("\nLoading GSE63061 (", round(file.info(f63061)$size/1e6, 1), "MB)...\n")

  gse63061 <- getGEO(filename = f63061, getGPL = FALSE)
  expr_63061 <- exprs(gse63061)
  meta_63061 <- pData(gse63061)

  cat("  Expression:", nrow(expr_63061), "probes x", ncol(expr_63061), "samples\n")

  ch1_cols_61 <- grep(":ch1$", colnames(meta_63061), value = TRUE)
  cat("  Metadata fields:\n")
  for (col in ch1_cols_61) {
    cat("    ", col, ":", paste(head(unique(meta_63061[[col]]), 5), collapse = " | "), "\n")
  }

  diag_col_61 <- NULL
  for (col in ch1_cols_61) {
    vals <- unique(tolower(meta_63061[[col]]))
    if (any(grepl("ad|alzheimer|control|mci|dementia", vals))) {
      diag_col_61 <- col
      cat("\n  Diagnosis column:", col, "\n")
      print(table(meta_63061[[col]]))
      break
    }
  }

  save(expr_63061, meta_63061, file = "results/GSE63061_processed.RData")
  cat("  Saved: results/GSE63061_processed.RData\n")
} else {
  cat("GSE63061 not available — skipping.\n")
}


###############################################################################
# SECTION 2: limma DEG on GSE63060 & GSE63061
#
# Both are Illumina HumanHT-12 BeadChip — same platform family as GSE140829.
# Processing: log2 → filter → quantile normalize → limma with covariates.
#
# IMPORTANT: Probe IDs are ILMN_* — need platform annotation for gene mapping.
# We'll download the GPL annotation separately or map via existing probe_gene_map.
###############################################################################

cat("\n\n=== SECTION 2: limma DEG on GSE63060/63061 ===\n\n")

run_limma_deg <- function(expr_mat, meta_df, diag_col, dataset_name,
                          age_col = NULL, sex_col = NULL) {
  # Determine groups
  diag_vals <- meta_df[[diag_col]]
  cat("  Raw diagnosis values:\n")
  print(table(diag_vals))

  # Standardize diagnosis labels
  diag_std <- tolower(trimws(diag_vals))
  diag_std <- gsub("alzheimer.*|alzheimer's disease", "ad", diag_std)
  diag_std <- gsub("^control$|^ctl$|^normal$|^healthy.*", "control", diag_std)
  diag_std <- gsub("^mci$|mild cognitive.*", "mci", diag_std)

  cat("  Standardized:\n")
  print(table(diag_std))

  # Keep only AD, Control, MCI
  keep_idx <- diag_std %in% c("ad", "control", "mci")
  if (sum(keep_idx) == 0) {
    cat("  WARNING: No AD/Control/MCI samples found. Trying original labels...\n")
    # Try with original case
    keep_idx <- tolower(diag_vals) %in% c("ad", "control", "mci",
                                           "alzheimer's disease", "ctl",
                                           "normal", "healthy")
    diag_std[tolower(diag_vals) %in% c("alzheimer's disease", "ad")] <- "ad"
    diag_std[tolower(diag_vals) %in% c("ctl", "control", "normal", "healthy")] <- "control"
    diag_std[tolower(diag_vals) %in% c("mci")] <- "mci"
  }

  expr_sub <- expr_mat[, keep_idx]
  meta_sub <- meta_df[keep_idx, ]
  meta_sub$group <- factor(diag_std[keep_idx], levels = c("control", "mci", "ad"))

  cat("  Samples for analysis:", ncol(expr_sub), "\n")
  cat("  Group sizes:\n")
  print(table(meta_sub$group))

  # Check if data needs log2 transform
  if (max(expr_sub[1:100, 1:min(5, ncol(expr_sub))], na.rm = TRUE) > 100) {
    cat("  Applying log2 transform...\n")
    expr_sub <- log2(pmax(expr_sub, 1))
  }

  # Filter low-expression probes (bottom 25%)
  med_expr <- apply(expr_sub, 1, median, na.rm = TRUE)
  keep_probes <- med_expr > quantile(med_expr, 0.25, na.rm = TRUE)
  expr_sub <- expr_sub[keep_probes, ]
  cat("  Probes after filtering:", nrow(expr_sub), "\n")

  # Quantile normalize
  expr_sub <- normalizeBetweenArrays(expr_sub, method = "quantile")

  # Build design matrix
  design <- model.matrix(~ 0 + group, data = meta_sub)
  colnames(design) <- levels(meta_sub$group)

  # Add covariates if available
  if (!is.null(age_col) && age_col %in% colnames(meta_sub)) {
    meta_sub$age_num <- as.numeric(meta_sub[[age_col]])
    if (sum(!is.na(meta_sub$age_num)) > ncol(expr_sub) * 0.5) {
      design <- model.matrix(~ 0 + group + age_num, data = meta_sub)
      colnames(design)[1:3] <- levels(meta_sub$group)
      cat("  Design: group + age\n")
    }
  }
  if (!is.null(sex_col) && sex_col %in% colnames(meta_sub)) {
    meta_sub$sex_bin <- as.numeric(factor(meta_sub[[sex_col]])) - 1
    if (length(unique(meta_sub$sex_bin)) == 2) {
      design <- cbind(design, sex_bin = meta_sub$sex_bin)
      cat("  Design: + sex\n")
    }
  }

  # Fit limma
  fit <- lmFit(expr_sub, design)

  contrast_matrix <- makeContrasts(
    AD_vs_Control = ad - control,
    MCI_vs_Control = mci - control,
    AD_vs_MCI = ad - mci,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # Extract AD vs Control
  deg <- topTable(fit2, coef = "AD_vs_Control", number = Inf, sort.by = "P")
  deg$probe_id <- rownames(deg)

  cat("\n  AD vs Control results:\n")
  cat("    Total probes:", nrow(deg), "\n")
  cat("    Significant (adj.P < 0.05):", sum(deg$adj.P.Val < 0.05), "\n")
  cat("    Up:", sum(deg$adj.P.Val < 0.05 & deg$logFC > 0), "\n")
  cat("    Down:", sum(deg$adj.P.Val < 0.05 & deg$logFC < 0), "\n")

  # Also get MCI vs Control
  deg_mci <- topTable(fit2, coef = "MCI_vs_Control", number = Inf, sort.by = "P")
  deg_mci$probe_id <- rownames(deg_mci)
  cat("  MCI vs Control significant:", sum(deg_mci$adj.P.Val < 0.05), "\n")

  return(list(
    AD_vs_Control = deg,
    MCI_vs_Control = deg_mci,
    n_ad = sum(meta_sub$group == "ad"),
    n_control = sum(meta_sub$group == "control"),
    n_mci = sum(meta_sub$group == "mci"),
    expr = expr_sub,
    meta = meta_sub
  ))
}

# --- Run DEG for GSE63060 ---
if (exists("expr_63060") && exists("meta_63060")) {
  cat("--- GSE63060 DEG Analysis ---\n")

  # Find the right column names
  age_col_60 <- grep("age", ch1_cols_60, value = TRUE, ignore.case = TRUE)[1]
  sex_col_60 <- grep("sex|gender", ch1_cols_60, value = TRUE, ignore.case = TRUE)[1]

  deg_63060 <- run_limma_deg(expr_63060, meta_63060, diag_col_60,
                              "GSE63060", age_col_60, sex_col_60)

  write.csv(deg_63060$AD_vs_Control,
            "results/DEG_GSE63060_AD_vs_Control.csv", row.names = FALSE)
  write.csv(deg_63060$MCI_vs_Control,
            "results/DEG_GSE63060_MCI_vs_Control.csv", row.names = FALSE)
  cat("  Saved DEG results.\n\n")
}

# --- Run DEG for GSE63061 ---
if (exists("expr_63061") && exists("meta_63061")) {
  cat("--- GSE63061 DEG Analysis ---\n")

  age_col_61 <- grep("age", ch1_cols_61, value = TRUE, ignore.case = TRUE)[1]
  sex_col_61 <- grep("sex|gender", ch1_cols_61, value = TRUE, ignore.case = TRUE)[1]

  deg_63061 <- run_limma_deg(expr_63061, meta_63061, diag_col_61,
                              "GSE63061", age_col_61, sex_col_61)

  write.csv(deg_63061$AD_vs_Control,
            "results/DEG_GSE63061_AD_vs_Control.csv", row.names = FALSE)
  write.csv(deg_63061$MCI_vs_Control,
            "results/DEG_GSE63061_MCI_vs_Control.csv", row.names = FALSE)
  cat("  Saved DEG results.\n\n")
}


###############################################################################
# SECTION 3: Map Probes to Genes & Standardize Effect Sizes
#
# GSE63060/63061 use ILMN_* probe IDs. We need to map to gene symbols.
# Strategy: download GPL annotation or use the probe_gene_map from GSE140829
# (same Illumina platform family — many probes overlap).
#
# For meta-analysis, we need per-gene:
#   - logFC (effect size)
#   - SE (standard error of logFC)
#   - study label
###############################################################################

cat("\n\n=== SECTION 3: Standardize Effect Sizes ===\n\n")

# Helper: collapse probe-level DEG results to gene level
# Uses the probe with the smallest p-value per gene
collapse_probes_to_genes <- function(deg_df, probe_to_gene) {
  # Merge with gene mapping
  merged <- merge(deg_df, probe_to_gene, by = "probe_id")
  merged <- merged[!is.na(merged$symbol) & merged$symbol != "", ]

  # Keep best probe per gene (lowest P.Value)
  merged <- merged[order(merged$symbol, merged$P.Value), ]
  merged <- merged[!duplicated(merged$symbol), ]
  rownames(merged) <- merged$symbol

  return(merged)
}

# Build a combined probe-to-gene mapping
# Start with GSE140829 mapping (already have it)
cat("Probe-gene mappings available from GSE140829:", nrow(probe_gene_map), "\n")

# For GSE63060/63061, try to get platform annotation
# These are typically GPL10558 (HumanHT-12 v4) or GPL6947 (HumanHT-12 v3)
# We can extract from the ExpressionSet annotation if available

get_probe_gene_map_from_eset <- function(eset, dataset_name) {
  # Try fData first
  fd <- fData(eset)
  if (ncol(fd) > 0) {
    cat("  ", dataset_name, "fData columns:", paste(head(colnames(fd), 10), collapse = ", "), "\n")

    # Look for gene symbol column
    sym_col <- grep("symbol|gene.symbol|gene_symbol|Symbol|SYMBOL",
                    colnames(fd), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(sym_col)) {
      map <- data.frame(
        probe_id = rownames(fd),
        symbol = as.character(fd[[sym_col]]),
        stringsAsFactors = FALSE
      )
      map <- map[!is.na(map$symbol) & map$symbol != "", ]
      cat("  Found", nrow(map), "probe-gene mappings from fData\n")
      return(map)
    }
  }

  # If fData is empty (because getGPL=FALSE), use GSE140829 mapping
  # ILMN probes are shared across Illumina platforms
  cat("  No fData — using GSE140829 probe-gene map (shared ILMN probes)\n")
  return(NULL)
}

# Try to get gene mappings for 63060/63061
map_63060 <- NULL
map_63061 <- NULL

if (exists("gse63060")) {
  map_63060 <- get_probe_gene_map_from_eset(gse63060, "GSE63060")
}
if (exists("gse63061")) {
  map_63061 <- get_probe_gene_map_from_eset(gse63061, "GSE63061")
}

# If we couldn't get maps from fData, download GPL annotation
if (is.null(map_63060) || is.null(map_63061)) {
  cat("\nDownloading GPL annotation for probe-to-gene mapping...\n")

  # Determine GPL from probe IDs
  if (exists("expr_63060")) {
    sample_probes <- head(rownames(expr_63060), 5)
    cat("Sample probe IDs from GSE63060:", paste(sample_probes, collapse = ", "), "\n")
  }

  # Try to download GPL10558 (most common Illumina HumanHT-12 v4)
  tryCatch({
    gpl <- getGEO("GPL10558", destdir = raw_dir)
    gpl_table <- Table(gpl)
    cat("GPL10558 table:", nrow(gpl_table), "rows\n")
    cat("Columns:", paste(head(colnames(gpl_table), 10), collapse = ", "), "\n")

    sym_col <- grep("symbol|Symbol|SYMBOL|gene_symbol",
                    colnames(gpl_table), value = TRUE, ignore.case = TRUE)[1]
    id_col <- grep("^ID$|^ID_REF$|^ILMN_ID$",
                   colnames(gpl_table), value = TRUE, ignore.case = TRUE)[1]

    if (!is.na(sym_col) && !is.na(id_col)) {
      gpl_map <- data.frame(
        probe_id = as.character(gpl_table[[id_col]]),
        symbol = as.character(gpl_table[[sym_col]]),
        stringsAsFactors = FALSE
      )
      gpl_map <- gpl_map[!is.na(gpl_map$symbol) & gpl_map$symbol != "", ]
      cat("GPL10558 mappings:", nrow(gpl_map), "\n")

      if (is.null(map_63060)) map_63060 <- gpl_map
      if (is.null(map_63061)) map_63061 <- gpl_map
    }
  }, error = function(e) {
    cat("GPL download failed:", e$message, "\n")
    cat("Falling back to GSE140829 probe map (shared ILMN probes).\n")
    # Use the GSE140829 map — it covers most ILMN probes
  })

  # Ultimate fallback: use GSE140829 map
  if (is.null(map_63060)) map_63060 <- probe_gene_map
  if (is.null(map_63061)) map_63061 <- probe_gene_map
}

# Now collapse all DEG results to gene level
cat("\n--- Collapsing to gene level ---\n")

# GSE140829: already at gene level (from Step 2)
cat("GSE140829: already gene-level,", nrow(deg_AD_ctrl), "genes\n")

# GSE63060
if (exists("deg_63060")) {
  deg_63060_genes <- collapse_probes_to_genes(deg_63060$AD_vs_Control, map_63060)
  cat("GSE63060:", nrow(deg_63060_genes), "genes after collapse\n")
}

# GSE63061
if (exists("deg_63061")) {
  deg_63061_genes <- collapse_probes_to_genes(deg_63061$AD_vs_Control, map_63061)
  cat("GSE63061:", nrow(deg_63061_genes), "genes after collapse\n")
}

# GSE270454: already at gene level
cat("GSE270454:", nrow(res_270_ad_nd), "genes\n")

# Build standardized effect size table
# Each row = one gene-study pair: gene, study, logFC, SE, P, N_ad, N_ctrl
cat("\n--- Building standardized effect size table ---\n")

effects_list <- list()

# Study 1: GSE140829 (limma)
# SE = logFC / t
se_140829 <- deg_AD_ctrl$logFC / deg_AD_ctrl$t
effects_list[["GSE140829"]] <- data.frame(
  gene = deg_AD_ctrl$gene,
  study = "GSE140829",
  platform = "microarray",
  population = "Japanese",
  logFC = deg_AD_ctrl$logFC,
  SE = se_140829,
  pvalue = deg_AD_ctrl$P.Value,
  adj_p = deg_AD_ctrl$adj.P.Val,
  n_case = 204,  # AD
  n_control = 249,
  stringsAsFactors = FALSE
)

# Study 2: GSE270454 (DESeq2) — use lfcSE directly
effects_list[["GSE270454"]] <- data.frame(
  gene = res_270_ad_nd$gene,
  study = "GSE270454",
  platform = "rnaseq",
  population = "Unknown",
  logFC = res_270_ad_nd$log2FoldChange,
  SE = res_270_ad_nd$lfcSE,
  pvalue = res_270_ad_nd$pvalue,
  adj_p = res_270_ad_nd$padj,
  n_case = 10,  # AD
  n_control = 25,  # ASM+ASO
  stringsAsFactors = FALSE
)

# Study 3: GSE63060 (limma)
if (exists("deg_63060_genes")) {
  se_63060 <- deg_63060_genes$logFC / deg_63060_genes$t
  effects_list[["GSE63060"]] <- data.frame(
    gene = deg_63060_genes$symbol,
    study = "GSE63060",
    platform = "microarray",
    population = "European",
    logFC = deg_63060_genes$logFC,
    SE = se_63060,
    pvalue = deg_63060_genes$P.Value,
    adj_p = deg_63060_genes$adj.P.Val,
    n_case = deg_63060$n_ad,
    n_control = deg_63060$n_control,
    stringsAsFactors = FALSE
  )
}

# Study 4: GSE63061 (limma)
if (exists("deg_63061_genes")) {
  se_63061 <- deg_63061_genes$logFC / deg_63061_genes$t
  effects_list[["GSE63061"]] <- data.frame(
    gene = deg_63061_genes$symbol,
    study = "GSE63061",
    platform = "microarray",
    population = "European",
    logFC = deg_63061_genes$logFC,
    SE = se_63061,
    pvalue = deg_63061_genes$P.Value,
    adj_p = deg_63061_genes$adj.P.Val,
    n_case = deg_63061$n_ad,
    n_control = deg_63061$n_control,
    stringsAsFactors = FALSE
  )
}

# Combine all
effects_all <- do.call(rbind, effects_list)

# Remove rows with NA/Inf SE
effects_all <- effects_all[is.finite(effects_all$SE) & effects_all$SE > 0 &
                            is.finite(effects_all$logFC) &
                            !is.na(effects_all$pvalue), ]

cat("\nCombined effect sizes:", nrow(effects_all), "gene-study pairs\n")
cat("Studies:", paste(unique(effects_all$study), collapse = ", "), "\n")
cat("Genes per study:\n")
print(table(effects_all$study))

# How many genes appear in 2+ studies?
gene_counts <- table(effects_all$gene)
cat("\nGenes in 1 study:", sum(gene_counts == 1), "\n")
cat("Genes in 2 studies:", sum(gene_counts == 2), "\n")
cat("Genes in 3 studies:", sum(gene_counts == 3), "\n")
cat("Genes in 4 studies:", sum(gene_counts == 4), "\n")

# Save
write.csv(effects_all, "results/standardized_effects.csv", row.names = FALSE)
cat("Saved: results/standardized_effects.csv\n")


###############################################################################
# SECTION 4: Per-Gene Random Effects Meta-Analysis (metafor REML)
#
# For each gene present in ≥2 studies, fit a random-effects model:
#   yi = logFC, vi = SE^2, method = "REML"
#
# Output per gene: pooled effect, pooled SE, pooled P, I² heterogeneity,
#                  tau², Q-test p-value, number of studies
###############################################################################

cat("\n\n=== SECTION 4: Per-Gene Meta-Analysis ===\n\n")

if (!requireNamespace("metafor", quietly = TRUE)) {
  cat("Installing metafor...\n")
  install.packages("metafor", repos = "https://cloud.r-project.org")
}
library(metafor)

# Only meta-analyze genes in ≥2 studies
genes_multi <- names(gene_counts[gene_counts >= 2])
cat("Genes in ≥2 studies:", length(genes_multi), "\n")
cat("Running per-gene REML meta-analysis...\n\n")

# Pre-allocate results
meta_results <- data.frame(
  gene = character(length(genes_multi)),
  n_studies = integer(length(genes_multi)),
  pooled_logFC = numeric(length(genes_multi)),
  pooled_SE = numeric(length(genes_multi)),
  pooled_pval = numeric(length(genes_multi)),
  I2 = numeric(length(genes_multi)),
  tau2 = numeric(length(genes_multi)),
  Q_pval = numeric(length(genes_multi)),
  direction_consistent = logical(length(genes_multi)),
  studies = character(length(genes_multi)),
  stringsAsFactors = FALSE
)

pb <- txtProgressBar(min = 0, max = length(genes_multi), style = 3)

for (i in seq_along(genes_multi)) {
  g <- genes_multi[i]
  dat <- effects_all[effects_all$gene == g, ]

  # Variance = SE^2
  vi <- dat$SE^2

  tryCatch({
    # Fit random-effects model
    fit <- rma(yi = dat$logFC, vi = vi, method = "REML")

    meta_results$gene[i] <- g
    meta_results$n_studies[i] <- nrow(dat)
    meta_results$pooled_logFC[i] <- fit$beta[1]
    meta_results$pooled_SE[i] <- fit$se
    meta_results$pooled_pval[i] <- fit$pval
    meta_results$I2[i] <- fit$I2
    meta_results$tau2[i] <- fit$tau2
    meta_results$Q_pval[i] <- fit$QEp
    meta_results$direction_consistent[i] <- length(unique(sign(dat$logFC))) == 1
    meta_results$studies[i] <- paste(dat$study, collapse = ";")

  }, error = function(e) {
    # Fallback: fixed-effect if REML fails (rare, e.g., 2 studies with same SE)
    tryCatch({
      fit <- rma(yi = dat$logFC, vi = vi, method = "FE")
      meta_results$gene[i] <<- g
      meta_results$n_studies[i] <<- nrow(dat)
      meta_results$pooled_logFC[i] <<- fit$beta[1]
      meta_results$pooled_SE[i] <<- fit$se
      meta_results$pooled_pval[i] <<- fit$pval
      meta_results$I2[i] <<- 0
      meta_results$tau2[i] <<- 0
      meta_results$Q_pval[i] <<- 1
      meta_results$direction_consistent[i] <<- length(unique(sign(dat$logFC))) == 1
      meta_results$studies[i] <<- paste(dat$study, collapse = ";")
    }, error = function(e2) {
      meta_results$gene[i] <<- g
      meta_results$n_studies[i] <<- nrow(dat)
    })
  })

  setTxtProgressBar(pb, i)
}
close(pb)

# Remove failed rows
meta_results <- meta_results[meta_results$gene != "" & !is.na(meta_results$pooled_pval), ]

# FDR correction
meta_results$adj_pval <- p.adjust(meta_results$pooled_pval, method = "BH")

# Sort by p-value
meta_results <- meta_results[order(meta_results$pooled_pval), ]

cat("\n\nMeta-analysis complete!\n")
cat("Total genes analyzed:", nrow(meta_results), "\n")
cat("Significant (adj P < 0.05):", sum(meta_results$adj_pval < 0.05, na.rm = TRUE), "\n")
cat("Direction consistent:", sum(meta_results$direction_consistent, na.rm = TRUE),
    "(", round(mean(meta_results$direction_consistent, na.rm = TRUE) * 100, 1), "%)\n")
cat("High heterogeneity (I² > 75%):", sum(meta_results$I2 > 75, na.rm = TRUE), "\n")

cat("\nTop 20 meta-analysis genes:\n")
print(head(meta_results[, c("gene", "n_studies", "pooled_logFC", "pooled_pval",
                             "adj_pval", "I2", "direction_consistent")], 20))

write.csv(meta_results, "results/meta_analysis_results.csv", row.names = FALSE)
cat("\nSaved: results/meta_analysis_results.csv\n")


###############################################################################
# SECTION 5: RS (Robustness Scoring)
#
# RS = PCA-driven composite across multiple robustness dimensions:
#
#   D1: Cross-Study Consistency  — direction concordance across studies
#   D2: Effect Size Magnitude    — |pooled logFC|
#   D3: Statistical Strength     — -log10(pooled P)
#   D4: Heterogeneity (inverse)  — (100 - I²)/100
#   D5: Cross-Platform Agreement — concordance between microarray & RNA-seq
#   D6: Replication Breadth      — n_studies / max_studies
#   D7: Stage Specificity Score  — MCI→AD progression signal
#
# Each dimension is scaled 0-1, then combined via PCA (first PC).
# Final RS is rescaled to [0, 100].
###############################################################################

cat("\n\n=== SECTION 5: Robustness Scoring (RS) ===\n\n")

# Start with meta-analysis results
rs_data <- meta_results[, c("gene", "n_studies", "pooled_logFC", "pooled_pval",
                             "adj_pval", "I2", "tau2", "direction_consistent")]

n_studies_max <- max(rs_data$n_studies, na.rm = TRUE)

# D1: Cross-Study Consistency (binary direction concordance → proportion)
# For each gene, what fraction of studies agree on direction?
cat("Computing D1: Cross-Study Consistency...\n")
rs_data$D1_consistency <- NA
for (i in seq_len(nrow(rs_data))) {
  g <- rs_data$gene[i]
  dat <- effects_all[effects_all$gene == g, ]
  if (nrow(dat) >= 2) {
    # Fraction of studies agreeing with the majority direction
    signs <- sign(dat$logFC)
    majority <- as.numeric(names(sort(table(signs), decreasing = TRUE))[1])
    rs_data$D1_consistency[i] <- mean(signs == majority)
  } else {
    rs_data$D1_consistency[i] <- 0.5  # uninformative for single-study
  }
}

# D2: Effect Size Magnitude (|pooled logFC|, capped at 99th percentile)
cat("Computing D2: Effect Size Magnitude...\n")
logfc_abs <- abs(rs_data$pooled_logFC)
cap <- quantile(logfc_abs, 0.99, na.rm = TRUE)
rs_data$D2_effect_size <- pmin(logfc_abs, cap) / cap

# D3: Statistical Strength (-log10 pooled P, capped)
cat("Computing D3: Statistical Strength...\n")
neglogp <- -log10(pmax(rs_data$pooled_pval, 1e-300))
cap_p <- quantile(neglogp, 0.99, na.rm = TRUE)
rs_data$D3_stat_strength <- pmin(neglogp, cap_p) / cap_p

# D4: Heterogeneity (inverse) — low I² = good
cat("Computing D4: Heterogeneity (inverse)...\n")
rs_data$D4_homogeneity <- (100 - pmin(rs_data$I2, 100)) / 100

# D5: Cross-Platform Agreement
# For genes present in both microarray AND RNA-seq studies,
# check direction concordance across platforms
cat("Computing D5: Cross-Platform Agreement...\n")
rs_data$D5_cross_platform <- NA
for (i in seq_len(nrow(rs_data))) {
  g <- rs_data$gene[i]
  dat <- effects_all[effects_all$gene == g, ]
  platforms <- unique(dat$platform)
  if (length(platforms) >= 2) {
    # Average sign per platform
    micro_sign <- sign(mean(dat$logFC[dat$platform == "microarray"]))
    rnaseq_sign <- sign(mean(dat$logFC[dat$platform == "rnaseq"]))
    rs_data$D5_cross_platform[i] <- ifelse(micro_sign == rnaseq_sign, 1.0, 0.0)
  } else {
    rs_data$D5_cross_platform[i] <- 0.5  # neutral if single platform
  }
}

# D6: Replication Breadth (n_studies / max)
cat("Computing D6: Replication Breadth...\n")
rs_data$D6_replication <- rs_data$n_studies / n_studies_max

# D7: Stage Specificity — MCI→AD progression signal
# Compare MCI vs Control effect with AD vs Control effect
# Genes showing AD > MCI signal = progressive biomarker
cat("Computing D7: Stage Specificity...\n")
rs_data$D7_stage_spec <- NA

# Get MCI vs Control results from GSE140829 (largest dataset)
for (i in seq_len(nrow(rs_data))) {
  g <- rs_data$gene[i]

  # AD vs Control logFC from meta
  ad_fc <- rs_data$pooled_logFC[i]

  # MCI vs Control logFC from GSE140829
  mci_idx <- which(deg_MCI_ctrl$gene == g)
  if (length(mci_idx) > 0) {
    mci_fc <- deg_MCI_ctrl$logFC[mci_idx[1]]

    # Stage specificity: same direction AND AD effect > MCI effect
    if (sign(ad_fc) == sign(mci_fc) && abs(ad_fc) > abs(mci_fc)) {
      # Progressive: score based on how much stronger AD signal is
      rs_data$D7_stage_spec[i] <- min(abs(ad_fc) / max(abs(mci_fc), 0.01), 3) / 3
    } else if (sign(ad_fc) == sign(mci_fc)) {
      # Same direction but MCI stronger — still somewhat informative
      rs_data$D7_stage_spec[i] <- 0.3
    } else {
      # Opposite directions — poor stage specificity
      rs_data$D7_stage_spec[i] <- 0.0
    }
  } else {
    rs_data$D7_stage_spec[i] <- 0.5  # no MCI data
  }
}

# Combine dimensions via PCA
cat("\nRunning PCA on robustness dimensions...\n")
dim_cols <- c("D1_consistency", "D2_effect_size", "D3_stat_strength",
              "D4_homogeneity", "D5_cross_platform", "D6_replication",
              "D7_stage_spec")

# Replace NAs with column medians for PCA
dim_matrix <- rs_data[, dim_cols]
for (col in dim_cols) {
  na_idx <- is.na(dim_matrix[[col]])
  if (any(na_idx)) {
    dim_matrix[[col]][na_idx] <- median(dim_matrix[[col]], na.rm = TRUE)
  }
}

# Scale each dimension 0-1 (some already are, but normalize for PCA)
for (col in dim_cols) {
  rng <- range(dim_matrix[[col]], na.rm = TRUE)
  if (rng[2] > rng[1]) {
    dim_matrix[[col]] <- (dim_matrix[[col]] - rng[1]) / (rng[2] - rng[1])
  }
}

# PCA
pca_fit <- prcomp(dim_matrix, center = TRUE, scale. = TRUE)
cat("PCA variance explained:\n")
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100
print(round(var_explained[1:min(5, length(var_explained))], 1))

cat("\nPC1 loadings (dimension weights):\n")
loadings <- pca_fit$rotation[, 1]
print(round(loadings, 3))

# RS = PC1 score, rescaled to [0, 100]
pc1_scores <- pca_fit$x[, 1]

# Flip if negatively correlated with D3 (we want higher = more robust)
if (cor(pc1_scores, dim_matrix$D3_stat_strength) < 0) {
  pc1_scores <- -pc1_scores
}

# Rescale to 0-100
rs_data$RS <- (pc1_scores - min(pc1_scores)) / (max(pc1_scores) - min(pc1_scores)) * 100

cat("\nRS distribution:\n")
cat("  Mean:", round(mean(rs_data$RS), 1), "\n")
cat("  Median:", round(median(rs_data$RS), 1), "\n")
cat("  SD:", round(sd(rs_data$RS), 1), "\n")
cat("  Range:", round(range(rs_data$RS), 1), "\n")

cat("\nTop 20 genes by RS:\n")
rs_top <- rs_data[order(-rs_data$RS), ]
print(head(rs_top[, c("gene", "RS", "n_studies", "pooled_logFC", "adj_pval",
                        "I2", "D1_consistency", "D5_cross_platform")], 20))

write.csv(rs_data, "results/RS_scores.csv", row.names = FALSE)
cat("\nSaved: results/RS_scores.csv\n")


###############################################################################
# SECTION 6: FM (Failure Mode Classification)
#
# Each gene is classified into its PRIMARY failure mode:
#
#   FM1 — Platform Dependency:    opposite direction across microarray vs RNA-seq
#   FM2 — Population Bias:        significant in one population but not others
#   FM3 — Stage Specificity:      significant in AD vs Ctrl but not MCI vs Ctrl
#                                 (or vice versa) — stage-dependent
#   FM4 — Temporal Instability:   high heterogeneity (I² > 75%) across studies
#   FM5 — Analytical Variability: (placeholder — requires DESeq2 vs edgeR vs
#                                  limma-voom comparison, computed separately)
#
# Genes with RS > 70 and no dominant failure mode → "Robust"
# Multiple FM flags can co-occur; we report the primary (strongest) one.
###############################################################################

cat("\n\n=== SECTION 6: Failure Mode Classification (FM) ===\n\n")

rs_data$FM1_platform <- FALSE
rs_data$FM2_population <- FALSE
rs_data$FM3_stage <- FALSE
rs_data$FM4_temporal <- FALSE
rs_data$FM5_analytical <- FALSE  # placeholder
rs_data$primary_FM <- "Robust"

for (i in seq_len(nrow(rs_data))) {
  g <- rs_data$gene[i]
  dat <- effects_all[effects_all$gene == g, ]

  # --- FM1: Platform Dependency ---
  platforms <- unique(dat$platform)
  if (length(platforms) >= 2) {
    micro_fc <- mean(dat$logFC[dat$platform == "microarray"], na.rm = TRUE)
    rnaseq_fc <- mean(dat$logFC[dat$platform == "rnaseq"], na.rm = TRUE)
    if (sign(micro_fc) != sign(rnaseq_fc) && abs(micro_fc) > 0.05 && abs(rnaseq_fc) > 0.05) {
      rs_data$FM1_platform[i] <- TRUE
    }
  }

  # --- FM2: Population Bias ---
  pops <- unique(dat$population)
  if (length(pops) >= 2) {
    # Check if significant in one population but not another
    pop_sig <- sapply(pops, function(p) {
      pop_dat <- dat[dat$population == p, ]
      any(pop_dat$pvalue < 0.05)
    })
    pop_nonsig <- sapply(pops, function(p) {
      pop_dat <- dat[dat$population == p, ]
      all(pop_dat$pvalue > 0.1)
    })
    # FM2 if significant in some populations but clearly non-significant in others
    if (any(pop_sig) && any(pop_nonsig)) {
      rs_data$FM2_population[i] <- TRUE
    }
  }

  # --- FM3: Stage Specificity ---
  # AD significant but MCI not (or vice versa) in GSE140829
  ad_idx <- which(deg_AD_ctrl$gene == g)
  mci_idx <- which(deg_MCI_ctrl$gene == g)
  if (length(ad_idx) > 0 && length(mci_idx) > 0) {
    ad_sig <- deg_AD_ctrl$adj.P.Val[ad_idx[1]] < 0.05
    mci_sig <- deg_MCI_ctrl$adj.P.Val[mci_idx[1]] < 0.05
    if (ad_sig != mci_sig) {
      rs_data$FM3_stage[i] <- TRUE
    }
  }

  # --- FM4: Temporal Instability (high I²) ---
  if (!is.na(rs_data$I2[i]) && rs_data$I2[i] > 75) {
    rs_data$FM4_temporal[i] <- TRUE
  }
}

# Assign primary FM (strongest signal)
for (i in seq_len(nrow(rs_data))) {
  fm_flags <- c(
    FM1 = rs_data$FM1_platform[i],
    FM2 = rs_data$FM2_population[i],
    FM3 = rs_data$FM3_stage[i],
    FM4 = rs_data$FM4_temporal[i]
  )

  if (any(fm_flags)) {
    # Priority: FM1 > FM4 > FM2 > FM3 (platform and temporal are most critical)
    if (fm_flags["FM1"]) {
      rs_data$primary_FM[i] <- "FM1_Platform"
    } else if (fm_flags["FM4"]) {
      rs_data$primary_FM[i] <- "FM4_Temporal"
    } else if (fm_flags["FM2"]) {
      rs_data$primary_FM[i] <- "FM2_Population"
    } else if (fm_flags["FM3"]) {
      rs_data$primary_FM[i] <- "FM3_Stage"
    }
  } else {
    rs_data$primary_FM[i] <- "Robust"
  }
}

cat("Failure Mode Distribution:\n")
print(table(rs_data$primary_FM))

cat("\nFM by RS quartile:\n")
rs_data$RS_quartile <- cut(rs_data$RS, breaks = quantile(rs_data$RS, c(0, 0.25, 0.5, 0.75, 1)),
                            labels = c("Q1_low", "Q2", "Q3", "Q4_high"), include.lowest = TRUE)
print(table(rs_data$RS_quartile, rs_data$primary_FM))

write.csv(rs_data, "results/RS_FM_scores.csv", row.names = FALSE)
cat("\nSaved: results/RS_FM_scores.csv\n")


###############################################################################
# SECTION 7: Summary Figures
###############################################################################

cat("\n\n=== SECTION 7: Figures ===\n\n")

# Fig 3: RS Distribution Histogram
cat("Generating Fig 3: RS Distribution...\n")
png("figures/Fig3_RS_Distribution.png", width = 800, height = 500)
par(mar = c(5, 5, 4, 2))
hist(rs_data$RS, breaks = 50, col = "steelblue", border = "white",
     main = "Robustness Score (RS) Distribution",
     xlab = "RS (0-100)", ylab = "Number of Genes",
     cex.lab = 1.3, cex.main = 1.4)
abline(v = median(rs_data$RS), lty = 2, lwd = 2, col = "red")
abline(v = 70, lty = 2, lwd = 2, col = "darkgreen")
legend("topright",
       legend = c(paste("Median:", round(median(rs_data$RS), 1)),
                  "Robust threshold (70)",
                  paste("N genes:", nrow(rs_data))),
       lty = c(2, 2, NA), lwd = c(2, 2, NA),
       col = c("red", "darkgreen", NA), bty = "n", cex = 1.1)
dev.off()
cat("Saved: figures/Fig3_RS_Distribution.png\n")

# Fig 4: FM Pie/Bar Chart
cat("Generating Fig 4: FM Distribution...\n")
fm_table <- table(rs_data$primary_FM)
fm_colors <- c("FM1_Platform" = "#e74c3c", "FM2_Population" = "#f39c12",
               "FM3_Stage" = "#3498db", "FM4_Temporal" = "#9b59b6",
               "Robust" = "#2ecc71")

png("figures/Fig4_FM_Distribution.png", width = 800, height = 500)
par(mar = c(6, 5, 4, 2))
bp <- barplot(fm_table[order(-fm_table)],
              col = fm_colors[names(fm_table[order(-fm_table)])],
              main = "Failure Mode Classification",
              ylab = "Number of Genes", las = 2,
              cex.lab = 1.3, cex.main = 1.4)
text(bp, fm_table[order(-fm_table)] + max(fm_table) * 0.02,
     labels = fm_table[order(-fm_table)], cex = 1.0)
dev.off()
cat("Saved: figures/Fig4_FM_Distribution.png\n")

# Fig 5: RS vs I² scatter (heterogeneity landscape)
cat("Generating Fig 5: RS vs I² Scatter...\n")
png("figures/Fig5_RS_vs_I2.png", width = 700, height = 700)
par(mar = c(5, 5, 4, 2))

fm_col <- fm_colors[rs_data$primary_FM]
fm_col[is.na(fm_col)] <- "grey50"

plot(rs_data$I2, rs_data$RS,
     pch = 20, cex = 0.6, col = adjustcolor(fm_col, alpha = 0.5),
     xlab = "I² Heterogeneity (%)", ylab = "Robustness Score (RS)",
     main = "RS vs Heterogeneity by Failure Mode",
     cex.lab = 1.3, cex.main = 1.4)
abline(h = 70, lty = 2, col = "darkgreen")
abline(v = 75, lty = 2, col = "red")

legend("topright",
       legend = names(fm_colors),
       col = fm_colors, pch = 20, cex = 0.9, pt.cex = 1.5,
       title = "Failure Mode")
dev.off()
cat("Saved: figures/Fig5_RS_vs_I2.png\n")

# Fig 6: Forest plot for top 10 genes
cat("Generating Fig 6: Forest plot (top 10 genes)...\n")
top10 <- head(rs_data[order(-rs_data$RS), ], 10)

png("figures/Fig6_Forest_Top10.png", width = 900, height = 700)
par(mar = c(5, 12, 4, 2))

n_genes <- nrow(top10)
y_pos <- n_genes:1

plot(NA, xlim = range(c(top10$pooled_logFC - 1.96 * meta_results$pooled_SE[match(top10$gene, meta_results$gene)],
                         top10$pooled_logFC + 1.96 * meta_results$pooled_SE[match(top10$gene, meta_results$gene)]),
                       na.rm = TRUE),
     ylim = c(0.5, n_genes + 0.5),
     xlab = "Pooled log2 Fold Change (AD vs Control)",
     ylab = "", yaxt = "n",
     main = "Top 10 Robust Biomarker Genes (by RS)")

abline(v = 0, lty = 2, col = "grey50")

for (j in seq_len(n_genes)) {
  gene_j <- top10$gene[j]
  fc <- top10$pooled_logFC[j]
  se <- meta_results$pooled_SE[match(gene_j, meta_results$gene)]

  # CI
  segments(fc - 1.96 * se, y_pos[j], fc + 1.96 * se, y_pos[j], lwd = 2)
  # Point
  points(fc, y_pos[j], pch = 18, cex = 1.5,
         col = ifelse(fc > 0, "red", "blue"))
}

axis(2, at = y_pos,
     labels = paste0(top10$gene, " (RS=", round(top10$RS, 0), ")"),
     las = 1, cex.axis = 0.9)

dev.off()
cat("Saved: figures/Fig6_Forest_Top10.png\n")


###############################################################################
# SECTION 8: Save Everything
###############################################################################

cat("\n\n=== Saving Full Workspace ===\n")
save.image("results/03_meta_RS_FM.RData")
cat("Saved: results/03_meta_RS_FM.RData\n")

cat("\n================================================================\n")
cat("  STEP 3 COMPLETE\n")
cat("================================================================\n\n")
cat("Key outputs:\n")
cat("  results/standardized_effects.csv  — all gene-study effect sizes\n")
cat("  results/meta_analysis_results.csv — per-gene pooled estimates\n")
cat("  results/RS_scores.csv             — robustness scores (7 dimensions)\n")
cat("  results/RS_FM_scores.csv          — RS + failure mode classification\n\n")
cat("  figures/Fig3_RS_Distribution.png   — RS histogram\n")
cat("  figures/Fig4_FM_Distribution.png   — FM bar chart\n")
cat("  figures/Fig5_RS_vs_I2.png          — RS vs heterogeneity\n")
cat("  figures/Fig6_Forest_Top10.png      — forest plot top 10\n\n")

cat("Robust genes (RS > 70):", sum(rs_data$RS > 70), "\n")
cat("High-confidence biomarkers (RS > 70, adj P < 0.05, direction consistent):\n")
hc <- rs_data[rs_data$RS > 70 & rs_data$adj_pval < 0.05 &
                rs_data$direction_consistent == TRUE, ]
hc <- hc[order(-hc$RS), ]
if (nrow(hc) > 0) {
  print(hc[, c("gene", "RS", "pooled_logFC", "adj_pval", "n_studies",
                "I2", "primary_FM")])
} else {
  cat("  (none at RS > 70 — try RS > 60)\n")
  hc <- rs_data[rs_data$RS > 60 & rs_data$adj_pval < 0.05 &
                  rs_data$direction_consistent == TRUE, ]
  hc <- hc[order(-hc$RS), ]
  print(head(hc[, c("gene", "RS", "pooled_logFC", "adj_pval", "n_studies",
                      "I2", "primary_FM")], 20))
}

cat("\n\nNext steps:\n")
cat("  1. Paper1_Step4_LOSO_CV.R — Leave-One-Study-Out cross-validation\n")
cat("  2. FM5 analytical variability (DESeq2 vs edgeR vs limma-voom)\n")
cat("  3. CIBERSORTx immune deconvolution\n")
cat("  4. RNA-to-Protein clinical translation mapping\n")
