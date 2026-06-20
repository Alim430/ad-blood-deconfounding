###############################################################################
# Step 13b (CORRECTED): ADNI explicit de-confounding validation.
# Mirrors Step13 EXACTLY (max-variance probe collapse, visit-matched DX,
# quantile-norm) and ONLY adds composition covariates, so the comparison is
# clean. Reports concordance of 162/143/356 sets WITHOUT vs WITH composition.
#   Rscript script/Paper1_Step13b_ADNI_deconfound.R
###############################################################################
suppressWarnings(suppressMessages({library(limma); library(ADNIMERGE)}))
setwd(".")
am <- get("adnimerge")

## ---- parse EXACTLY like Step13 ----------------------------------------------
raw <- read.csv("data/ADNI/ADNI_Gene_Expression_Profile.csv", header=FALSE,
                stringsAsFactors=FALSE, check.names=FALSE)
hdr <- which(raw[[1]]=="ProbeSet")[1]; samp_cols <- 4:ncol(raw)
visit <- tolower(as.character(raw[2, samp_cols])); subj <- as.character(raw[3, samp_cols])
sym  <- raw[(hdr+1):nrow(raw), 3]
expr <- as.matrix(raw[(hdr+1):nrow(raw), samp_cols]); mode(expr) <- "numeric"
rownames(expr) <- sym; colnames(expr) <- subj
expr <- expr[!is.na(sym) & sym!="" & sym!="---", ]
# collapse probes -> genes (max-variance probe), exactly as Step13
v <- apply(expr,1,var,na.rm=TRUE); o <- order(rownames(expr), -v)
expr <- expr[o,]; expr <- expr[!duplicated(rownames(expr)),]
cat(sprintf("ADNI: %d genes x %d samples\n", nrow(expr), ncol(expr)))

## ---- visit-matched DX (the correct assignment) ------------------------------
RID <- as.integer(sub(".*_S_","",subj))
meta <- data.frame(col=seq_along(subj), RID=RID, visit=visit, stringsAsFactors=FALSE)
dx_at <- function(rid,vc){ r<-am[am$RID==rid & am$VISCODE==vc,]; if(nrow(r)) r$DX[1] else NA }
dx_bl <- function(rid){ r<-am[am$RID==rid & am$VISCODE=="bl",]; if(nrow(r)) r$DX[1] else NA }
get1  <- function(rid,col){ r<-am[am$RID==rid,]; r<-r[!is.na(r[[col]]),]; if(nrow(r)) r[[col]][1] else NA }
meta$DX  <- mapply(function(r,v){ d<-dx_at(r,v); if(is.na(d)) d<-dx_bl(r); d }, meta$RID, meta$visit)
meta$AGE <- sapply(meta$RID, get1, "AGE")
meta$SEX <- sapply(meta$RID, function(r) as.character(get1(r,"PTGENDER")))
meta$APOE4 <- sapply(meta$RID, function(r) suppressWarnings(as.numeric(get1(r,"APOE4"))))
meta$group <- ifelse(meta$DX %in% c("Dementia","AD"),"AD",
              ifelse(meta$DX=="CN","Control", ifelse(meta$DX=="MCI","MCI",NA)))
keep <- !is.na(meta$group) & meta$group %in% c("Control","AD")
md <- meta[keep,]; X <- expr[, md$col, drop=FALSE]
X <- normalizeBetweenArrays(X, method="quantile")
md$grp <- factor(md$group, levels=c("Control","AD"))
md$AGE[is.na(md$AGE)] <- median(md$AGE, na.rm=TRUE); md$APOE4[is.na(md$APOE4)] <- 0
cat(sprintf("CN vs AD usable: %d (AD=%d, CN=%d)  [should match Step13: 343/93/250]\n",
            nrow(md), sum(md$grp=="AD"), sum(md$grp=="Control")))

## ---- composition on the SAME quantile-normed log2 matrix (no double log) -----
DANAHER <- list(
  Neutrophil=c("CSF3R","S100A12","FCAR","SIGLEC5","TNFRSF10C","FPR1","MMP25"),
  Monocyte=c("CD68","CD163","CD14","CTSS","FCN1","MS4A6A","TYROBP"),
  Bcell=c("CD19","CD79A","CD79B","MS4A1","TCL1A","TNFRSF13B"),
  Tcell_CD4=c("CD3D","CD3E","CD3G","CD4","TRAT1","ITM2A"),
  Tcell_CD8=c("CD8A","CD8B","GZMK"),
  NK=c("KLRD1","KLRK1","NCR1","NCAM1","FGFBP2","GNLY","NKG7"))
zsig <- function(m, sets){
  mz <- (m - rowMeans(m,na.rm=TRUE)) / pmax(apply(m,1,sd,na.rm=TRUE),1e-6)
  sapply(sets, function(g){ g<-intersect(g,rownames(mz))
    if(!length(g)) return(rep(NA_real_,ncol(mz))); colMeans(mz[g,,drop=FALSE],na.rm=TRUE) }) }
comp <- zsig(X, DANAHER)                          # samples x lineages, aligned to X cols = md rows
for (j in seq_len(ncol(comp))) comp[is.na(comp[,j]),j] <- median(comp[,j],na.rm=TRUE)
md <- cbind(md, comp)
cat("composition vs AD collinearity: ",
    paste(colnames(comp), sprintf("%+.2f",sapply(colnames(comp), function(c) cor(md[[c]], as.integer(md$grp=="AD")))),
          sep="=", collapse="  "), "\n")

## ---- two models: NO-comp (reproduce Step13) vs WITH-comp ---------------------
des0 <- model.matrix(~ grp + AGE + factor(SEX) + APOE4, data=md)
des1 <- model.matrix(~ grp + AGE + factor(SEX) + APOE4 +
                       Neutrophil+Monocyte+Bcell+Tcell_CD4+Tcell_CD8+NK, data=md)
tt0 <- topTable(eBayes(lmFit(X,des0)), coef="grpAD", number=Inf, sort.by="none"); tt0$gene<-rownames(X)
tt1 <- topTable(eBayes(lmFit(X,des1)), coef="grpAD", number=Inf, sort.by="none"); tt1$gene<-rownames(X)

## ---- concordance of survivor sets vs discovery direction --------------------
s <- read.csv("results/biomarker_survival_after_adjustment.csv")
disc_dir <- setNames(sign(s$logFC_adj), s$gene)
rdset <- function(f) intersect(readLines(f), rownames(X))
sets <- list("162_cohort_specific"=rdset("results/sensitivity_162.txt"),
             "143_robust_core"=rdset("results/robust_core_143.txt"),
             "356_uniform"=rdset("results/sensitivity_356.txt"))
conc <- function(ttx,g){ sub<-ttx[match(g,ttx$gene),]; k<-sum(sign(sub$logFC)==disc_dir[g],na.rm=TRUE)
  n<-length(g); list(n=n,k=k,pct=100*k/n,z=(k-0.5*n)/(0.5*sqrt(n)),rep=sum(sign(sub$logFC)==disc_dir[g] & sub$P.Value<0.05,na.rm=TRUE)) }
out <- data.frame()
cat("\n=== ADNI concordance with discovery direction ===\n")
for (nm in names(sets)){ g<-sets[[nm]]; a<-conc(tt0,g); b<-conc(tt1,g)
  out <- rbind(out, data.frame(set=nm, n=a$n,
     nocomp_pct=round(a$pct), nocomp_z=round(a$z,1), nocomp_rep=a$rep,
     withcomp_pct=round(b$pct), withcomp_z=round(b$z,1), withcomp_rep=b$rep))
  cat(sprintf("  %-20s n=%3d | NO-comp %2.0f%% (z=%+.1f, rep=%d) | WITH-comp %2.0f%% (z=%+.1f, rep=%d)\n",
              nm, a$n, a$pct, a$z, a$rep, b$pct, b$z, b$rep)) }
write.csv(out, "results/ADNI_deconfounded_concordance.csv", row.names=FALSE)
cat("\nKey test: does WITH-comp concordance stay HIGH (z>3)? If yes -> survivors replicate beyond composition.\n")
cat("Saved: results/ADNI_deconfounded_concordance.csv\n")
