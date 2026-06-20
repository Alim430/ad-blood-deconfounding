###############################################################################
# Step 6c: UNIFORM de-confounding sensitivity analysis.
# Re-run the adjustment with the SAME model in all three cohorts
# (~ group + age + sex + composition), i.e. only covariates available in ALL
# cohorts, then re-meta and compare survivors to the cohort-specific result.
# Tests whether using non-identical per-cohort models changed the conclusions.
#   Rscript script/Paper1_Step6c_uniform_deconfound.R
###############################################################################
suppressWarnings(suppressMessages({library(limma); library(metafor)}))
e <- new.env(); load("results/06_adjusted_reanalysis.RData", envir=e)
run_adj_limma <- get("run_adj_limma", e); expr_g <- get("expr_g", e)
cell_scores   <- get("cell_scores", e)
meta140 <- get("meta140",e); meta60 <- get("meta60",e); meta61 <- get("meta61",e)

# UNIFORM model: NO extra_terms (no batch/APOE) for any cohort -> composition+age+sex only
cat("Running UNIFORM adjusted DE (composition+age+sex) in all 3 cohorts...\n")
d140 <- run_adj_limma(expr_g$GSE140829, meta140, cell_scores$GSE140829[colnames(expr_g$GSE140829),], cohort="GSE140829")
d60  <- run_adj_limma(expr_g$GSE63060,  meta60,  cell_scores$GSE63060[colnames(expr_g$GSE63060),],  cohort="GSE63060")
d61  <- run_adj_limma(expr_g$GSE63061,  meta61,  cell_scores$GSE63061[colnames(expr_g$GSE63061),],  cohort="GSE63061")

eff <- rbind(data.frame(gene=d140$gene, logFC=d140$logFC, SE=d140$SE, study="GSE140829"),
             data.frame(gene=d60$gene,  logFC=d60$logFC,  SE=d60$SE,  study="GSE63060"),
             data.frame(gene=d61$gene,  logFC=d61$logFC,  SE=d61$SE,  study="GSE63061"))
eff <- eff[is.finite(eff$logFC) & is.finite(eff$SE) & eff$SE>0, ]
tab <- table(eff$gene); multi <- names(tab[tab>=2])
cat("Genes in >=2 cohorts:", length(multi), "\n  meta-analysing...\n")

res <- data.frame(gene=multi, pooled_logFC=NA_real_, pooled_pval=NA_real_, dir_consistent=NA)
for (i in seq_along(multi)){
  dat <- eff[eff$gene==multi[i],]
  fit <- tryCatch(rma(yi=dat$logFC, vi=dat$SE^2, method="REML"),
                  error=function(z) tryCatch(rma(yi=dat$logFC, vi=dat$SE^2, method="FE"), error=function(z2) NULL))
  if (is.null(fit)) next
  res$pooled_logFC[i] <- fit$beta[1]; res$pooled_pval[i] <- fit$pval
  res$dir_consistent[i] <- length(unique(sign(dat$logFC)))==1
}
res$adjP <- p.adjust(res$pooled_pval, "BH")

# compare to current survivors
s <- read.csv("results/biomarker_survival_after_adjustment.csv")
t1 <- s[s$tier_orig=="Tier1_Robust", c("gene","logFC_orig","survives")]
m <- merge(t1, res, by="gene")
m$survives_uniform <- !is.na(m$adjP) & m$adjP<0.05 & m$dir_consistent & sign(m$pooled_logFC)==sign(m$logFC_orig)
cur <- m$survives %in% c(TRUE,"TRUE","True")
cat(sprintf("\n=== Tier-1 (723) survival ===\n"))
cat(sprintf("  cohort-specific model (current): %d survive\n", sum(cur)))
cat(sprintf("  UNIFORM model (comp+age+sex)    : %d survive\n", sum(m$survives_uniform)))
cat(sprintf("  overlap (survive in BOTH)       : %d\n", sum(cur & m$survives_uniform)))
cat(sprintf("  Jaccard(current, uniform)       : %.2f\n",
            sum(cur & m$survives_uniform)/sum(cur | m$survives_uniform)))
write.csv(m, "results/uniform_deconfound_compare.csv", row.names=FALSE)
cat("\nSaved: results/uniform_deconfound_compare.csv\n")
