###############################################################################
# Paper 1 — Step 13: Independent validation in ADNI blood (Affymetrix U219)
#
# THE credibility step: does the de-confounded survivor signature replicate, and
# classify CN vs AD, in a large INDEPENDENT cohort (~745 ADNI blood samples)?
#
# INPUTS (already downloaded to data/ADNI/):
#   ADNI_Gene_Expression_Profile.csv   (quirky format: 8 metadata rows on top)
#   ADNIMERGE_0.0.1.tar.gz             (clinical: DX, AGE, PTGENDER, APOE4)
#   results/biomarker_survival_after_adjustment.csv (survivors + discovery dir)
#
# OUTPUTS:
#   results/ADNI_survivor_replication.csv  — per-gene replication in ADNI
#   results/ADNI_validation_summary.txt
#   figures/pub/Fig12_ADNI_score.pdf, Fig12_ADNI_ROC.pdf
###############################################################################
setwd(".")
suppressWarnings(suppressMessages({ library(limma) }))

## ---- clinical (ADNIMERGE) ----
if (!requireNamespace("ADNIMERGE", quietly=TRUE))
  install.packages("data/ADNI/ADNIMERGE_0.0.1.tar.gz", repos=NULL, type="source")
suppressWarnings(suppressMessages(library(ADNIMERGE)))
am <- get("adnimerge")
am$RID <- as.integer(as.character(am$RID))
am$VISCODE <- tolower(trimws(as.character(am$VISCODE)))
am$DX <- as.character(am$DX)

## ---- parse the ADNI gene-expression matrix ----
cat("Parsing ADNI expression (large, ~1 min)...\n")
raw <- read.csv("data/ADNI/ADNI_Gene_Expression_Profile.csv", header=FALSE,
                stringsAsFactors=FALSE, check.names=FALSE)
hdr <- which(raw[[1]]=="ProbeSet")[1]                  # row with column header
ncol_all <- ncol(raw)
# samples are columns 4..end; drop trailing empty cols
samp_cols <- 4:ncol_all
phase <- as.character(raw[1, samp_cols]); visit <- tolower(as.character(raw[2, samp_cols]))
subj  <- as.character(raw[3, samp_cols])
keepc <- !is.na(subj) & subj!=""
samp_cols <- samp_cols[keepc]; phase<-phase[keepc]; visit<-visit[keepc]; subj<-subj[keepc]
sym  <- raw[(hdr+1):nrow(raw), 3]
expr <- as.matrix(raw[(hdr+1):nrow(raw), samp_cols])
suppressWarnings(class(expr)<-"numeric"); expr <- matrix(as.numeric(expr), nrow=length(sym))
rownames(expr)<-sym; colnames(expr)<-subj
expr <- expr[!is.na(sym) & sym!="" & sym!="---", ]
cat(sprintf("  %d probes x %d samples\n", nrow(expr), ncol(expr)))

# collapse probes -> genes (max-variance probe)
v <- apply(expr,1,var,na.rm=TRUE)
o <- order(rownames(expr), -v); expr<-expr[o,]; expr<-expr[!duplicated(rownames(expr)),]
cat(sprintf("  %d genes after collapse\n", nrow(expr)))

## ---- sample metadata + diagnosis join ----
RID <- as.integer(sub(".*_S_","",subj))
meta <- data.frame(col=seq_along(subj), subj=subj, RID=RID, visit=visit, phase=phase,
                   stringsAsFactors=FALSE)
# DX at the matching visit; fallback to baseline (bl) for that RID
dx_at <- function(rid,vc){ r<-am[am$RID==rid & am$VISCODE==vc,]; if(nrow(r)) r$DX[1] else NA }
dx_bl <- function(rid){ r<-am[am$RID==rid & am$VISCODE=="bl",]; if(nrow(r)) r$DX[1] else NA }
get1 <- function(rid,col){ r<-am[am$RID==rid,]; r<-r[!is.na(r[[col]]),]; if(nrow(r)) r[[col]][1] else NA }
meta$DX <- mapply(function(r,v){ d<-dx_at(r,v); if(is.na(d)) d<-dx_bl(r); d}, meta$RID, meta$visit)
meta$AGE <- sapply(meta$RID, get1, "AGE")
meta$SEX <- sapply(meta$RID, function(r) as.character(get1(r,"PTGENDER")))
meta$APOE4 <- sapply(meta$RID, function(r) suppressWarnings(as.numeric(get1(r,"APOE4"))))
meta$group <- ifelse(meta$DX %in% c("Dementia","AD"),"AD",
              ifelse(meta$DX=="CN","Control", ifelse(meta$DX=="MCI","MCI",NA)))
cat("ADNI diagnosis distribution:\n"); print(table(meta$group, useNA="ifany"))

## ---- restrict to CN/AD with expression ----
keep <- !is.na(meta$group) & meta$group %in% c("Control","AD")
md <- meta[keep,]; X <- expr[, md$col, drop=FALSE]
X <- normalizeBetweenArrays(X, method="quantile")
md$grp <- factor(md$group, levels=c("Control","AD"))
cat(sprintf("\nCN vs AD usable: %d (AD=%d, CN=%d)\n", nrow(md), sum(md$grp=="AD"), sum(md$grp=="Control")))

## ---- survivors + discovery direction ----
s <- read.csv("results/biomarker_survival_after_adjustment.csv", stringsAsFactors=FALSE)
surv <- s[s$tier_orig=="Tier1_Robust" & s$survives %in% c("TRUE","True"),
          c("gene","logFC_adj")]
surv <- surv[surv$gene %in% rownames(X), ]
cat(sprintf("Survivors measurable in ADNI: %d / %d\n", nrow(surv),
            sum(s$tier_orig=="Tier1_Robust" & s$survives %in% c("TRUE","True"))))

## ---- (1) per-gene replication: limma AD vs CN in ADNI (adjust age/sex/APOE) ----
des <- model.matrix(~ grp + AGE + factor(SEX) + APOE4,
                    data=transform(md, AGE=ifelse(is.na(AGE),median(AGE,na.rm=TRUE),AGE),
                                   APOE4=ifelse(is.na(APOE4),0,APOE4)))
fit <- eBayes(lmFit(X, des))
tt <- topTable(fit, coef="grpAD", number=Inf, sort.by="none"); tt$gene<-rownames(X)
rep <- merge(surv, tt[,c("gene","logFC","P.Value","adj.P.Val")], by="gene")
rep$concordant <- sign(rep$logFC)==sign(rep$logFC_adj)
rep$replicated <- rep$concordant & rep$P.Value<0.05
write.csv(rep, "results/ADNI_survivor_replication.csv", row.names=FALSE)
cat(sprintf("\nReplication: %d/%d concordant direction (%.0f%%); %d also nominal p<0.05\n",
            sum(rep$concordant), nrow(rep), 100*mean(rep$concordant), sum(rep$replicated)))

## ---- (2) survivor SIGNATURE SCORE -> AUC (CN vs AD) ----
Z <- t(scale(t(X[surv$gene,,drop=FALSE]))); Z[is.na(Z)]<-0
score <- colSums(Z * sign(surv$logFC_adj)) / nrow(surv)      # sign-weighted
md$score <- score
if(!requireNamespace("pROC",quietly=TRUE)) install.packages("pROC",repos="https://cloud.r-project.org")
library(pROC)
roc1 <- roc(md$grp, md$score, quiet=TRUE)
cat(sprintf("Survivor-signature score AUC (CN vs AD): %.3f\n", as.numeric(auc(roc1))))

## ---- (3) elastic-net classifier on survivor panel (10-fold CV AUC) ----
if(!requireNamespace("glmnet",quietly=TRUE)) install.packages("glmnet",repos="https://cloud.r-project.org")
library(glmnet)
Xp <- t(scale(t(X[surv$gene,]))); Xp[is.na(Xp)]<-0
cv <- cv.glmnet(t(Xp), md$grp, family="binomial", alpha=.5, nfolds=10, type.measure="auc")
cat(sprintf("Elastic-net survivor panel CV AUC: %.3f\n", max(cv$cvm)))

## ---- figures ----
suppressWarnings(suppressMessages(library(ggplot2)))
th<-theme_classic(base_size=10)+theme(plot.title=element_text(face="bold"))
gcol<-c(Control="#0072B2",AD="#D55E00")
ggsave("figures/pub/Fig12_ADNI_score.pdf",
  ggplot(md,aes(grp,score,fill=grp))+geom_boxplot(outlier.size=.5)+geom_jitter(width=.15,size=.3,alpha=.3)+
    scale_fill_manual(values=gcol,guide="none")+labs(title=sprintf("ADNI: survivor signature (AUC=%.2f)",as.numeric(auc(roc1))),x=NULL,y="signature score")+th,
  width=70,height=75,units="mm")
rdf<-data.frame(spec=rev(roc1$specificities),sens=rev(roc1$sensitivities))
ggsave("figures/pub/Fig12_ADNI_ROC.pdf",
  ggplot(rdf,aes(1-spec,sens))+geom_abline(linetype=2,color="grey")+geom_line(color="#D55E00",linewidth=.7)+
    coord_equal()+labs(title=sprintf("ADNI CN vs AD (AUC=%.2f)",as.numeric(auc(roc1))),x="1 - specificity",y="sensitivity")+th,
  width=70,height=70,units="mm")

sink("results/ADNI_validation_summary.txt")
cat("ADNI independent validation\n");cat("CN vs AD n:",nrow(md),"\n")
cat("Survivors measurable:",nrow(surv),"\n")
cat(sprintf("Direction-concordant: %d/%d (%.0f%%); nominal-replicated: %d\n",sum(rep$concordant),nrow(rep),100*mean(rep$concordant),sum(rep$replicated)))
cat(sprintf("Signature-score AUC: %.3f\n",as.numeric(auc(roc1))))
cat(sprintf("Elastic-net panel CV AUC: %.3f\n",max(cv$cvm)))
sink()
cat("\nSaved: results/ADNI_survivor_replication.csv, ADNI_validation_summary.txt, figures/pub/Fig12_*\n")
