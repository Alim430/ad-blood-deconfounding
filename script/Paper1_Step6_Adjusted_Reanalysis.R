###############################################################################
# Paper 1 — Step 6: Confounder-Adjusted Re-Analysis  (supersedes Steps 2-4)
#
# WHY THIS EXISTS
#   The original pipeline (Steps 2-4) adjusted DEGs for age + sex only. Review
#   of the data showed three unused / mishandled confounders that can manufacture
#   "robust" biomarkers:
#     (1) BATCH      — GSE140829 has 19 batches, never modeled
#     (2) APOE       — GSE140829 has APOE genotype on 575/587, never used
#     (3) CELL COMP  — whole-blood cell-proportion shifts (neutrophil up,
#                      lymphocyte down in AD) were never adjusted; the top
#                      "biomarkers" (CD14, FPR1, SIGLEC5, CD79B...) are cell
#                      lineage markers, so the signal may be composition, not
#                      disease-intrinsic transcription.
#   Plus: GSE270454 (n=45, 1 DEG) was a co-equal meta member; here it is DEMOTED
#   to a qualitative cross-platform check only.
#
# WHAT THIS SCRIPT DOES
#   Section 1: Estimate blood cell-type scores per sample (marker-based,
#              Danaher et al. 2017 style — no external reference download).
#   Section 2: Re-run limma DEG with the ADJUSTED model:
#                GSE140829: ~ group + age + sex + batch + APOE_e4 + cell scores
#                GSE63060 : ~ group + age + sex + cell scores
#                GSE63061 : ~ group + age + sex + cell scores
#   Section 3: Re-run random-effects meta-analysis (metafor REML) on the
#              3 MICROARRAY cohorts only.
#   Section 4: KEY RESULT — biomarker survival: how many original Tier-1 genes
#              remain significant + direction-consistent after adjustment.
#   Section 5: Honest validation — leave-one-cohort-out (LOCO) elastic-net
#              classifier, AD vs Control, report out-of-cohort AUC.
#   Section 6: GSE270454 cross-platform concordance (demoted / supportive only).
#
# OUTPUTS (all suffixed _adj so nothing overwrites the originals):
#   results/cell_scores_<cohort>.csv
#   results/DEG_<cohort>_AD_vs_Control_adj.csv
#   results/meta_analysis_results_adj.csv
#   results/biomarker_survival_after_adjustment.csv     <-- the headline table
#   results/LOCO_AUC.csv
#   results/06_adjusted_reanalysis.RData
#
# PREREQUISITE: results/04_LOSO_FM5_final.RData   (cumulative workspace)
# Run in RStudio section by section.
###############################################################################

setwd(".")

library(limma)
library(data.table)

dir.create("results", showWarnings = FALSE)

cat("Loading cumulative workspace (results/04_LOSO_FM5_final.RData)...\n")
load("results/04_LOSO_FM5_final.RData")
cat("Loaded.\n\n")


###############################################################################
# HELPERS
###############################################################################

# Collapse a probe-level expression matrix to gene level using a probe->symbol
# map, keeping the most variable probe per gene (IQR).
collapse_expr_to_genes <- function(expr_mat, probe2sym) {
  pmap <- probe2sym[probe2sym$probe_id %in% rownames(expr_mat), ]
  pmap <- pmap[!is.na(pmap$symbol) & pmap$symbol != "", ]
  iqrs <- apply(expr_mat[pmap$probe_id, , drop = FALSE], 1, IQR, na.rm = TRUE)
  ord  <- order(pmap$symbol, -iqrs)
  pmap <- pmap[ord, ]
  pmap <- pmap[!duplicated(pmap$symbol), ]
  out  <- expr_mat[pmap$probe_id, , drop = FALSE]
  rownames(out) <- pmap$symbol
  out
}

# Danaher-style cell-type score: mean log2 expression of marker genes present,
# z-scored across samples so each lineage contributes on a common scale.
# Markers chosen to be lineage-specific in whole blood.
CELL_MARKERS <- list(
  Neutrophil  = c("FCGR3B","CSF3R","CXCR2","CEACAM3","VNN3","FFAR2","MME","CSF2RB"),
  Monocyte    = c("CD14","CSF1R","FCN1","LYZ","CD163","MS4A7","SERPINA1","VCAN"),
  Bcell       = c("CD19","CD79A","MS4A1","CD22","TCL1A","BANK1","FCRL1"),
  Tcell_CD4   = c("CD3D","CD3E","CD3G","IL7R","CD4","CD28","CD5"),
  Tcell_CD8   = c("CD8A","CD8B","GZMK","CD2"),
  NK          = c("KLRD1","NKG7","GNLY","KLRF1","NCR1","NCAM1","KLRB1")
)

cell_scores_from_expr <- function(expr_genes_mat) {
  # expr_genes_mat: genes (symbols) x samples, log2 scale
  z <- t(scale(t(expr_genes_mat)))           # z-score each gene across samples
  z[is.na(z)] <- 0
  sc <- sapply(names(CELL_MARKERS), function(ct) {
    gs <- intersect(CELL_MARKERS[[ct]], rownames(z))
    if (length(gs) == 0) return(rep(NA_real_, ncol(z)))
    colMeans(z[gs, , drop = FALSE], na.rm = TRUE)
  })
  sc <- as.data.frame(sc)
  rownames(sc) <- colnames(expr_genes_mat)
  attr(sc, "n_markers") <- sapply(names(CELL_MARKERS),
                                  function(ct) length(intersect(CELL_MARKERS[[ct]], rownames(z))))
  sc
}

# Build gene-level expression for each cohort from whatever is in the workspace.
get_gene_expr <- function(cohort) {
  if (cohort == "GSE140829") {
    if (exists("expr_genes")) return(expr_genes)
    stopifnot(exists("expr_140829_matched"), exists("probe_gene_map"))
    e <- expr_140829_matched
    if (max(e[1:50, 1:5], na.rm = TRUE) > 100) e <- log2(e + 1)
    return(collapse_expr_to_genes(e, probe_gene_map))
  }
  if (cohort == "GSE63060") {
    stopifnot(exists("expr_63060"), exists("map_63060"))
    e <- expr_63060
    if (max(e[1:50, 1:5], na.rm = TRUE) > 100) e <- log2(pmax(e, 1))
    return(collapse_expr_to_genes(e, map_63060))
  }
  if (cohort == "GSE63061") {
    stopifnot(exists("expr_63061"), exists("map_63061"))
    e <- expr_63061
    if (max(e[1:50, 1:5], na.rm = TRUE) > 100) e <- log2(pmax(e, 1))
    return(collapse_expr_to_genes(e, map_63061))
  }
  stop("unknown cohort")
}


###############################################################################
# SECTION 1: Cell-type composition scores per cohort
###############################################################################

cat("=== SECTION 1: Cell-composition scores ===\n\n")

expr_g <- list(
  GSE140829 = get_gene_expr("GSE140829"),
  GSE63060  = get_gene_expr("GSE63060"),
  GSE63061  = get_gene_expr("GSE63061")
)

cell_scores <- lapply(names(expr_g), function(co) {
  cs <- cell_scores_from_expr(expr_g[[co]])
  cat(co, "- markers found per lineage:\n"); print(attr(cs, "n_markers"))
  write.csv(cbind(sample = rownames(cs), cs),
            sprintf("results/cell_scores_%s.csv", co), row.names = FALSE)
  cs
})
names(cell_scores) <- names(expr_g)
cat("\nSaved per-cohort cell scores.\n\n")


###############################################################################
# SECTION 2: Adjusted limma DEG
#
# Returns gene-level AD-vs-Control topTable with moderated SE (= logFC / t).
###############################################################################

cat("=== SECTION 2: Adjusted DEG ===\n\n")

# --- assemble per-cohort meta with group + covariates ---

# GSE140829 -----------------------------------------------------------------
meta140 <- meta_140829_matched
meta140$group <- factor(meta140$diagnosis, levels = c("Control","MCI","AD"))
meta140$age   <- suppressWarnings(as.numeric(meta140$age_at_draw))
meta140$sex   <- factor(meta140$Sex)
meta140$batch <- factor(meta140$batch)
# APOE: count of E4 alleles (0/1/2); missing -> NA then median-impute (=0) + flag
e4 <- sapply(strsplit(as.character(meta140$apoe_status), "_"),
             function(x) sum(x == "E4"))
e4[meta140$apoe_status %in% c("", NA)] <- NA
meta140$apoe_e4 <- e4
meta140$apoe_miss <- as.integer(is.na(meta140$apoe_e4))
meta140$apoe_e4[is.na(meta140$apoe_e4)] <- 0

# Robustly reorder meta140 so its rows match the expression columns, trying
# every plausible key (geo_accession, beadchip_id, rownames, title-embedded id).
align_meta_to_expr <- function(expr_mat, meta_df) {
  cn <- colnames(expr_mat)
  keys <- list()
  for (k in c("geo_accession","beadchip_id","title")) {
    if (k %in% colnames(meta_df)) keys[[k]] <- as.character(meta_df[[k]])
  }
  keys[["rownames"]] <- rownames(meta_df)
  for (kn in names(keys)) {
    kv <- keys[[kn]]
    if (sum(cn %in% kv) == length(cn) && !any(duplicated(kv[kv %in% cn]))) {
      cat("  Aligned GSE140829 meta to expr via:", kn, "\n")
      return(meta_df[match(cn, kv), , drop = FALSE])
    }
  }
  # last resort: same N and already in order
  if (nrow(meta_df) == length(cn)) {
    cat("  WARNING: no key matched; assuming expr cols and meta rows are in the SAME ORDER.\n")
    return(meta_df)
  }
  stop(sprintf(paste0("Cannot align expr (%d cols) to meta (%d rows).\n",
                      "  expr cols: %s\n  geo: %s\n  beadchip: %s"),
               length(cn), nrow(meta_df),
               paste(head(cn,3),collapse=", "),
               paste(head(meta_df$geo_accession,3),collapse=", "),
               paste(head(meta_df$beadchip_id,3),collapse=", ")))
}
meta140 <- align_meta_to_expr(expr_g$GSE140829, meta140)
rownames(meta140) <- colnames(expr_g$GSE140829)
cs140 <- cell_scores$GSE140829[colnames(expr_g$GSE140829), ]

# GSE63060 / GSE63061 -------------------------------------------------------
# pull the :ch1 fields parsed in Step 3 (meta_63060 / meta_63061)
prep_anm <- function(meta_df, expr_mat) {
  # Only search the sample-characteristic fields (":ch1"), NOT generic GEO
  # columns like "status" (= "Public on ..."), which would shadow "status:ch1".
  ch1cols <- grep(":ch1$", colnames(meta_df), value = TRUE)
  ch <- function(p) {
    col <- grep(p, ch1cols, value = TRUE, ignore.case = TRUE)[1]
    if (is.na(col)) return(rep(NA, nrow(meta_df)))
    meta_df[[col]]
  }
  status <- tolower(trimws(ch("status")))
  grp <- ifelse(status %in% c("ad","alzheimer's disease"), "AD",
         ifelse(status %in% c("ctl","control","normal"), "Control",
         ifelse(status == "mci", "MCI", NA)))
  data.frame(
    group = factor(grp, levels = c("Control","MCI","AD")),
    age   = suppressWarnings(as.numeric(ch("age"))),
    sex   = factor(ch("gender|sex")),
    row.names = rownames(meta_df),
    stringsAsFactors = FALSE
  )
}
meta60 <- prep_anm(meta_63060, expr_g$GSE63060)
meta61 <- prep_anm(meta_63061, expr_g$GSE63061)

# Run adjusted limma for one cohort and return AD-vs-Control gene table.
run_adj_limma <- function(expr_genes_mat, meta_df, cell_sc, extra_terms = NULL,
                          cohort = "") {
  keep <- !is.na(meta_df$group) & meta_df$group %in% c("Control","MCI","AD")
  e  <- expr_genes_mat[, keep, drop = FALSE]
  md <- droplevels(meta_df[keep, , drop = FALSE])
  cs <- cell_sc[keep, , drop = FALSE]

  # quantile normalize gene matrix
  e <- normalizeBetweenArrays(e, method = "quantile")

  # base design: 0 + group so we can contrast AD - Control
  md$group <- factor(md$group, levels = c("Control","MCI","AD"))
  df <- data.frame(group = md$group, age = md$age, sex = md$sex)
  # cell scores
  df <- cbind(df, cs)
  # cohort-specific extra covariates (batch, apoe) supplied via extra_terms df
  if (!is.null(extra_terms)) df <- cbind(df, extra_terms[keep, , drop = FALSE])

  # drop covariate columns that are constant or fully NA
  covs <- setdiff(colnames(df), "group")
  usable <- covs[sapply(covs, function(c) {
    v <- df[[c]]
    if (is.factor(v)) return(nlevels(droplevels(v)) > 1)
    length(unique(v[!is.na(v)])) > 1 && mean(is.na(v)) < 0.5
  })]
  # median-impute numeric NAs so no samples are dropped
  for (c in usable) if (is.numeric(df[[c]])) df[[c]][is.na(df[[c]])] <- median(df[[c]], na.rm = TRUE)

  form <- as.formula(paste("~ 0 + group +", paste(usable, collapse = " + ")))
  design <- model.matrix(form, data = df)
  # rename the three group columns to plain names for contrasts
  colnames(design)[grep("^groupControl$", colnames(design))] <- "Control"
  colnames(design)[grep("^groupMCI$",     colnames(design))] <- "MCI"
  colnames(design)[grep("^groupAD$",      colnames(design))] <- "AD"

  cat(sprintf("  %s: %d samples | covariates: %s\n",
              cohort, ncol(e), paste(usable, collapse = ", ")))

  fit  <- lmFit(e, design)
  cm   <- makeContrasts(AD_vs_Control = AD - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tt   <- topTable(fit2, coef = "AD_vs_Control", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt$SE   <- tt$logFC / tt$t            # moderated SE
  tt
}

cat("Running adjusted DEG (this includes batch + APOE + cell scores)...\n")
deg140_adj <- run_adj_limma(expr_g$GSE140829, meta140, cs140,
                            extra_terms = meta140[, c("batch","apoe_e4","apoe_miss")],
                            cohort = "GSE140829")
deg60_adj  <- run_adj_limma(expr_g$GSE63060, meta60,
                            cell_scores$GSE63060[colnames(expr_g$GSE63060), ],
                            cohort = "GSE63060")
deg61_adj  <- run_adj_limma(expr_g$GSE63061, meta61,
                            cell_scores$GSE63061[colnames(expr_g$GSE63061), ],
                            cohort = "GSE63061")

for (nm in c("deg140_adj","deg60_adj","deg61_adj")) {
  d <- get(nm)
  co <- c(deg140_adj="GSE140829", deg60_adj="GSE63060", deg61_adj="GSE63061")[nm]
  cat(sprintf("  %s adjusted: %d sig (adj.P<0.05)\n", co, sum(d$adj.P.Val < 0.05, na.rm = TRUE)))
  write.csv(d, sprintf("results/DEG_%s_AD_vs_Control_adj.csv", co), row.names = FALSE)
}
cat("\n")


###############################################################################
# SECTION 3: Adjusted meta-analysis (3 microarray cohorts only)
###############################################################################

cat("=== SECTION 3: Adjusted meta-analysis (microarray only) ===\n\n")

if (!requireNamespace("metafor", quietly = TRUE))
  install.packages("metafor", repos = "https://cloud.r-project.org")
library(metafor)

mk_eff <- function(d, study) data.frame(
  gene = d$gene, study = study, logFC = d$logFC, SE = d$SE,
  pvalue = d$P.Value, stringsAsFactors = FALSE)

eff_adj <- rbind(mk_eff(deg140_adj, "GSE140829"),
                 mk_eff(deg60_adj,  "GSE63060"),
                 mk_eff(deg61_adj,  "GSE63061"))
eff_adj <- eff_adj[is.finite(eff_adj$SE) & eff_adj$SE > 0 &
                   is.finite(eff_adj$logFC) & !is.na(eff_adj$pvalue), ]

gene_n <- table(eff_adj$gene)
genes_multi <- names(gene_n[gene_n >= 2])
cat("Genes in >=2 microarray cohorts:", length(genes_multi), "\n")

meta_adj <- data.frame(gene = genes_multi, n_studies = NA_integer_,
                       pooled_logFC = NA_real_, pooled_SE = NA_real_,
                       pooled_pval = NA_real_, I2 = NA_real_,
                       direction_consistent = NA, stringsAsFactors = FALSE)
pb <- txtProgressBar(0, length(genes_multi), style = 3)
for (i in seq_along(genes_multi)) {
  dat <- eff_adj[eff_adj$gene == genes_multi[i], ]
  ok <- tryCatch({
    fit <- rma(yi = dat$logFC, vi = dat$SE^2, method = "REML"); TRUE
  }, error = function(e) {
    fit <<- tryCatch(rma(yi = dat$logFC, vi = dat$SE^2, method = "FE"),
                     error = function(e2) NULL); !is.null(fit)
  })
  if (ok) {
    meta_adj$n_studies[i]    <- nrow(dat)
    meta_adj$pooled_logFC[i] <- fit$beta[1]
    meta_adj$pooled_SE[i]    <- fit$se
    meta_adj$pooled_pval[i]  <- fit$pval
    meta_adj$I2[i]           <- ifelse(is.null(fit$I2), 0, fit$I2)
    meta_adj$direction_consistent[i] <- length(unique(sign(dat$logFC))) == 1
  }
  setTxtProgressBar(pb, i)
}
close(pb)
meta_adj <- meta_adj[!is.na(meta_adj$pooled_pval), ]
meta_adj$adj_pval <- p.adjust(meta_adj$pooled_pval, method = "BH")
meta_adj <- meta_adj[order(meta_adj$pooled_pval), ]

cat("\nAdjusted meta: ", nrow(meta_adj), "genes | sig(adj.P<0.05):",
    sum(meta_adj$adj_pval < 0.05), "\n")
write.csv(meta_adj, "results/meta_analysis_results_adj.csv", row.names = FALSE)


###############################################################################
# SECTION 4: HEADLINE — biomarker survival after adjustment
###############################################################################

cat("\n=== SECTION 4: Biomarker survival after adjustment ===\n\n")

# original final table (unadjusted) is `final_table` from Step 4
orig <- final_table[, c("gene","RS","tier","pooled_logFC","adj_pval","primary_FM")]
colnames(orig) <- c("gene","RS_orig","tier_orig","logFC_orig","adjP_orig","FM_orig")

surv <- merge(orig, meta_adj[, c("gene","pooled_logFC","adj_pval","I2",
                                 "direction_consistent")],
              by = "gene", all.x = TRUE)
colnames(surv)[colnames(surv)=="pooled_logFC"] <- "logFC_adj"
colnames(surv)[colnames(surv)=="adj_pval"]     <- "adjP_adj"

surv$survives <- !is.na(surv$adjP_adj) & surv$adjP_adj < 0.05 &
                 sign(surv$logFC_adj) == sign(surv$logFC_orig)
surv$sign_flip <- !is.na(surv$logFC_adj) &
                  sign(surv$logFC_adj) != sign(surv$logFC_orig)

t1 <- surv[surv$tier_orig == "Tier1_Robust", ]
cat(sprintf("Original Tier-1 genes: %d\n", nrow(t1)))
cat(sprintf("  Still in adjusted meta: %d\n", sum(!is.na(t1$adjP_adj))))
cat(sprintf("  SURVIVE (sig + same direction): %d (%.0f%%)\n",
            sum(t1$survives), 100*mean(t1$survives)))
cat(sprintf("  Direction FLIPPED after adjustment: %d\n", sum(t1$sign_flip, na.rm = TRUE)))

cat("\nTop original Tier-1 genes and their fate:\n")
show <- t1[order(-t1$RS_orig),
           c("gene","RS_orig","FM_orig","logFC_orig","adjP_orig",
             "logFC_adj","adjP_adj","survives")]
print(head(show, 30), row.names = FALSE)

write.csv(surv[order(-surv$RS_orig), ],
          "results/biomarker_survival_after_adjustment.csv", row.names = FALSE)
cat("\nSaved: results/biomarker_survival_after_adjustment.csv\n")


###############################################################################
# SECTION 5: Honest validation — Leave-One-Cohort-Out classifier AUC
#
# Train elastic-net logistic (AD vs Control) on 2 microarray cohorts,
# test AUC on the held-out cohort. Rotate over all 3. This is a REAL
# out-of-cohort generalization estimate (replaces the pseudo-LOSO).
###############################################################################

cat("\n=== SECTION 5: Leave-One-Cohort-Out AUC ===\n\n")

if (!requireNamespace("glmnet", quietly = TRUE))
  install.packages("glmnet", repos = "https://cloud.r-project.org")
if (!requireNamespace("pROC", quietly = TRUE))
  install.packages("pROC", repos = "https://cloud.r-project.org")
library(glmnet); library(pROC)

# Build AD-vs-Control gene x sample matrices + labels per cohort
build_xy <- function(expr_genes_mat, grp) {
  keep <- grp %in% c("Control","AD")
  list(X = t(expr_genes_mat[, keep, drop = FALSE]),
       y = factor(ifelse(grp[keep] == "AD", 1, 0)))
}
xy <- list(
  GSE140829 = build_xy(expr_g$GSE140829, as.character(meta140$group)),
  GSE63060  = build_xy(expr_g$GSE63060,  as.character(meta60$group)),
  GSE63061  = build_xy(expr_g$GSE63061,  as.character(meta61$group))
)

# restrict to genes common to all three cohorts, and to top-200 adjusted-meta genes
common_genes <- Reduce(intersect, lapply(xy, function(z) colnames(z$X)))
panel <- head(meta_adj$gene[meta_adj$gene %in% common_genes], 200)
cat("Classifier panel size:", length(panel), "genes\n\n")

loco <- data.frame(held_out = names(xy), auc = NA_real_, n_test = NA_integer_)
for (k in seq_along(xy)) {
  test_co  <- names(xy)[k]
  train_co <- setdiff(names(xy), test_co)
  Xtr <- do.call(rbind, lapply(train_co, function(c) scale(xy[[c]]$X[, panel, drop = FALSE])))
  ytr <- unlist(lapply(train_co, function(c) xy[[c]]$y))
  Xte <- scale(xy[[test_co]]$X[, panel, drop = FALSE])
  yte <- xy[[test_co]]$y
  Xtr[is.na(Xtr)] <- 0; Xte[is.na(Xte)] <- 0
  cvf <- cv.glmnet(as.matrix(Xtr), ytr, family = "binomial", alpha = 0.5)
  pr  <- as.numeric(predict(cvf, as.matrix(Xte), s = "lambda.min", type = "response"))
  ro  <- roc(yte, pr, quiet = TRUE)
  loco$auc[k]    <- as.numeric(auc(ro))
  loco$n_test[k] <- length(yte)
  cat(sprintf("  Hold out %s: AUC = %.3f (n=%d)\n", test_co, loco$auc[k], loco$n_test[k]))
}
cat(sprintf("\nMean LOCO AUC: %.3f\n", mean(loco$auc)))
write.csv(loco, "results/LOCO_AUC.csv", row.names = FALSE)


###############################################################################
# SECTION 6: GSE270454 cross-platform check (DEMOTED — supportive only)
###############################################################################

cat("\n=== SECTION 6: GSE270454 cross-platform concordance (supportive) ===\n\n")

if (exists("res_270_ad_nd")) {
  cp <- merge(meta_adj[, c("gene","pooled_logFC")],
              res_270_ad_nd[, c("gene","log2FoldChange")], by = "gene")
  cp <- cp[is.finite(cp$pooled_logFC) & is.finite(cp$log2FoldChange), ]
  r  <- cor(cp$pooled_logFC, cp$log2FoldChange, use = "complete.obs")
  conc <- mean(sign(cp$pooled_logFC) == sign(cp$log2FoldChange)) * 100
  cat(sprintf("Genes compared: %d | r = %.3f | direction concordance = %.1f%%\n",
              nrow(cp), r, conc))
  cat("NOTE: GSE270454 is n=45 (10 AD) with 1 DEG — interpret as qualitative only.\n")
} else cat("res_270_ad_nd not in workspace — skipping.\n")


###############################################################################
# SAVE
###############################################################################

save.image("results/06_adjusted_reanalysis.RData")
cat("\n================================================================\n")
cat("  STEP 6 COMPLETE — adjusted re-analysis saved\n")
cat("================================================================\n")
cat("Headline: results/biomarker_survival_after_adjustment.csv\n")
cat("Validation: results/LOCO_AUC.csv (real out-of-cohort AUC)\n")
cat("Saved workspace: results/06_adjusted_reanalysis.RData\n")
