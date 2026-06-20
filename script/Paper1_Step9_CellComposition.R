###############################################################################
# Paper 1 — Step 9: Blood Cell-Composition Shift across Control/MCI/AD [LAYER 4]
#
# PURPOSE
#   The headline of the reframed paper is that most "robust" AD blood biomarkers
#   are CELL-COMPOSITION-MEDIATED. This step quantifies that composition shift
#   directly: do neutrophil/monocyte scores rise and lymphocyte (B/T/NK) scores
#   fall along Control -> MCI -> AD? This is the mechanism, and it validates why
#   adjusting for cell scores (Step 6) removed 78% of the signal.
#
#   Deconvolution-lite: reuses the 6 marker-based cell scores already computed in
#   Step 6 (Danaher-style), so NO new download/package. Tests each lineage across
#   diagnosis with Kruskal-Wallis (3-group), Spearman trend (CN<MCI<AD ordinal),
#   and AD-vs-Control Wilcoxon, per cohort, then summarizes direction across the
#   three microarray cohorts.
#   (Upgrade path for the final figure: CIBERSORTx / MCPcounter for true % — but
#    the marker-score test is sufficient and defensible for the claim.)
#
# PREREQUISITE: results/06_adjusted_reanalysis.RData (has cell_scores, meta140/60/61)
# OUTPUTS:
#   results/cell_composition_shift.csv   — per cohort x lineage stats
#   results/cell_composition_summary.csv — cross-cohort direction summary
#   figures/Fig_cell_composition_shift.png
###############################################################################

setwd(".")
dir.create("figures", showWarnings = FALSE)

load("results/06_adjusted_reanalysis.RData")
cat("Loaded Step-6 workspace.\n")

lineages <- c("Neutrophil","Monocyte","Bcell","Tcell_CD4","Tcell_CD8","NK")
myeloid  <- c("Neutrophil","Monocyte")
lymphoid <- c("Bcell","Tcell_CD4","Tcell_CD8","NK")

# group vector per cohort (already aligned to cell_scores row order in Step 6)
groups <- list(GSE140829 = meta140$group, GSE63060 = meta60$group, GSE63061 = meta61$group)

# assemble long table: cohort, sample, group, lineage, score
rows <- list()
for (co in names(cell_scores)) {
  cs <- cell_scores[[co]]
  g  <- groups[[co]]
  stopifnot(nrow(cs) == length(g))
  keep <- !is.na(g) & g %in% c("Control","MCI","AD")
  cs <- cs[keep, , drop = FALSE]; g <- droplevels(factor(g[keep], levels = c("Control","MCI","AD")))
  for (ln in lineages) {
    rows[[paste(co, ln)]] <- data.frame(cohort = co, group = g,
                                        lineage = ln, score = cs[[ln]],
                                        stringsAsFactors = FALSE)
  }
}
long <- do.call(rbind, rows)
long$group <- factor(long$group, levels = c("Control","MCI","AD"))

# ---- per cohort x lineage statistics ----
stat_rows <- list()
for (co in unique(long$cohort)) for (ln in lineages) {
  d <- long[long$cohort == co & long$lineage == ln, ]
  ord <- as.integer(d$group) - 1L                       # Control=0, MCI=1, AD=2
  kw  <- tryCatch(kruskal.test(score ~ group, data = d)$p.value, error = function(e) NA)
  sp  <- suppressWarnings(cor.test(d$score, ord, method = "spearman"))
  ad  <- d$score[d$group == "AD"]; cn <- d$score[d$group == "Control"]
  wt  <- tryCatch(wilcox.test(ad, cn)$p.value, error = function(e) NA)
  stat_rows[[paste(co, ln)]] <- data.frame(
    cohort = co, lineage = ln,
    n = nrow(d),
    KW_p = kw,
    trend_rho = unname(sp$estimate), trend_p = sp$p.value,
    AD_minus_CN = mean(ad, na.rm = TRUE) - mean(cn, na.rm = TRUE),
    AD_vs_CN_p = wt, stringsAsFactors = FALSE)
}
stats <- do.call(rbind, stat_rows)
stats <- stats[order(stats$lineage, stats$cohort), ]
write.csv(stats, "results/cell_composition_shift.csv", row.names = FALSE)

cat("\n=== Per-cohort cell-composition shift (AD vs Control) ===\n")
print(stats[, c("cohort","lineage","AD_minus_CN","AD_vs_CN_p","trend_rho","trend_p")],
      row.names = FALSE, digits = 3)

# ---- cross-cohort summary: is the shift consistent? ----
summ_rows <- list()
for (ln in lineages) {
  s <- stats[stats$lineage == ln, ]
  summ_rows[[ln]] <- data.frame(
    lineage = ln,
    mean_AD_minus_CN = mean(s$AD_minus_CN),
    n_cohorts_up   = sum(s$AD_minus_CN > 0),
    n_cohorts_down = sum(s$AD_minus_CN < 0),
    consistent = (all(s$AD_minus_CN > 0) | all(s$AD_minus_CN < 0)),
    n_sig_trend = sum(s$trend_p < 0.05, na.rm = TRUE),
    class = ifelse(ln %in% myeloid, "myeloid", "lymphoid"),
    stringsAsFactors = FALSE)
}
summary_tab <- do.call(rbind, summ_rows)
summary_tab <- summary_tab[order(-summary_tab$mean_AD_minus_CN), ]
write.csv(summary_tab, "results/cell_composition_summary.csv", row.names = FALSE)

cat("\n=== Cross-cohort summary (does the composition story hold?) ===\n")
print(summary_tab, row.names = FALSE, digits = 3)

up_myeloid   <- summary_tab[summary_tab$class=="myeloid"   & summary_tab$mean_AD_minus_CN > 0, "lineage"]
dn_lymphoid  <- summary_tab[summary_tab$class=="lymphoid"  & summary_tab$mean_AD_minus_CN < 0, "lineage"]
cat(sprintf("\nMyeloid UP in AD: %s\nLymphoid DOWN in AD: %s\n",
            paste(up_myeloid, collapse=", "), paste(dn_lymphoid, collapse=", ")))
cat("If myeloid up + lymphoid down across cohorts, the composition-mediated thesis holds.\n")

# ---- figure: 6 lineage panels, score by group (pooled across cohorts; scores are
#      within-cohort z-scored so comparable) ----
png("figures/Fig_cell_composition_shift.png", width = 1100, height = 720, res = 120)
par(mfrow = c(2, 3), mar = c(3.5, 4, 3, 1))
gcols <- c(Control = "#4C72B0", MCI = "#DD8452", AD = "#C44E52")
for (ln in lineages) {
  d <- long[long$lineage == ln, ]
  boxplot(score ~ group, data = d, col = gcols, outline = FALSE,
          main = ln, xlab = "", ylab = "cell-score (z)", border = "grey30")
  # trend annotation (pooled spearman)
  ord <- as.integer(d$group) - 1L
  rho <- suppressWarnings(cor(d$score, ord, method = "spearman", use = "complete.obs"))
  mtext(sprintf("trend rho = %+.2f", rho), side = 3, line = -1.2, cex = 0.7,
        col = ifelse(abs(rho) > 0.1, "black", "grey60"))
}
dev.off()
cat("\nSaved: figures/Fig_cell_composition_shift.png\n")
cat("Saved: results/cell_composition_shift.csv, results/cell_composition_summary.csv\n")
