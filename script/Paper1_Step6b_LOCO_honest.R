###############################################################################
# Step 6b: HONEST leave-one-cohort-out AUC (fixes feature-selection leakage).
# Original LOCO chose the 200-gene panel from the ALL-cohort adjusted meta,
# including the held-out cohort -> optimistic. Here, within each fold the panel
# is re-selected from ONLY the two training cohorts' adjusted DE (Stouffer Z).
# Also reports an all-common-genes elastic-net (no preselection) as a 2nd honest
# estimate. Loads the saved Step6 workspace; no GEO re-download.
#   Rscript script/Paper1_Step6b_LOCO_honest.R
# Output: results/LOCO_AUC_honest.csv
###############################################################################
suppressWarnings(suppressMessages({library(glmnet); library(pROC)}))
setwd(".")
e <- new.env(); load("results/06_adjusted_reanalysis.RData", envir=e)
expr_g <- get("expr_g",e)
grp <- list(GSE140829=as.character(get("meta140",e)$group),
            GSE63060 =as.character(get("meta60", e)$group),
            GSE63061 =as.character(get("meta61", e)$group))
deg <- list(GSE140829=get("deg140_adj",e), GSE63060=get("deg60_adj",e), GSE63061=get("deg61_adj",e))

build_xy <- function(m,g){ keep<-g %in% c("Control","AD")
  list(X=t(m[,keep,drop=FALSE]), y=factor(ifelse(g[keep]=="AD",1,0))) }
xy <- Map(build_xy, expr_g, grp)
common <- Reduce(intersect, lapply(xy, function(z) colnames(z$X)))
cat(sprintf("Common genes across 3 cohorts: %d\n\n", length(common)))

# training-only feature ranking: Stouffer Z over the 2 training cohorts' adjusted DE
rank_train <- function(train_co){
  zs <- lapply(train_co, function(c){
    d <- deg[[c]]; d <- d[d$gene %in% common, ]
    p <- pmax(pmin(d$P.Value,1),1e-300)
    setNames(sign(d$logFC)*qnorm(1-p/2), d$gene) })
  g <- Reduce(intersect, lapply(zs, names))
  Z <- Reduce(`+`, lapply(zs, function(z) z[g]))/sqrt(length(train_co))
  names(sort(abs(Z), decreasing=TRUE)) }

fit_auc <- function(panel, train_co, test_co){
  Xtr <- do.call(rbind, lapply(train_co, function(c) scale(xy[[c]]$X[,panel,drop=FALSE])))
  ytr <- unlist(lapply(train_co, function(c) xy[[c]]$y))
  Xte <- scale(xy[[test_co]]$X[,panel,drop=FALSE]); yte <- xy[[test_co]]$y
  Xtr[is.na(Xtr)]<-0; Xte[is.na(Xte)]<-0
  cvf <- cv.glmnet(as.matrix(Xtr), ytr, family="binomial", alpha=0.5)
  pr  <- as.numeric(predict(cvf, as.matrix(Xte), s="lambda.min", type="response"))
  as.numeric(auc(roc(yte, pr, quiet=TRUE))) }

out <- data.frame()
for (test_co in names(xy)){
  train_co <- setdiff(names(xy), test_co)
  panel200 <- head(rank_train(train_co), 200)          # leak-free top-200
  a_sel <- fit_auc(panel200, train_co, test_co)
  a_all <- fit_auc(common,   train_co, test_co)         # no preselection
  out <- rbind(out, data.frame(held_out=test_co, n_test=length(xy[[test_co]]$y),
                               auc_trainSelected200=round(a_sel,3),
                               auc_allGenes=round(a_all,3)))
  cat(sprintf("  Hold out %-10s  top200(train-only) AUC=%.3f   allGenes AUC=%.3f  (n=%d)\n",
              test_co, a_sel, a_all, length(xy[[test_co]]$y)))
}
out <- rbind(out, data.frame(held_out="MEAN", n_test=NA,
   auc_trainSelected200=round(mean(out$auc_trainSelected200),3),
   auc_allGenes=round(mean(out$auc_allGenes),3)))
write.csv(out, "results/LOCO_AUC_honest.csv", row.names=FALSE)
cat("\n(leaky original was 0.69 / 0.84 / 0.83)\nSaved: results/LOCO_AUC_honest.csv\n")
print(out, row.names=FALSE)
