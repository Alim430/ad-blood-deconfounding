###############################################################################
# AI Module 1a: export the survivor-panel expression matrix from ADNI for ML
# Output: results/ai/adni_panel.csv  (rows=samples; cols=genes + group/age/sex/apoe)
###############################################################################
setwd(".")
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(ADNIMERGE)))
am <- get("adnimerge"); am$RID<-as.integer(as.character(am$RID))
am$VISCODE<-tolower(trimws(as.character(am$VISCODE))); am$DX<-as.character(am$DX)

raw <- read.csv("data/ADNI/ADNI_Gene_Expression_Profile.csv", header=FALSE,
                stringsAsFactors=FALSE, check.names=FALSE)
hdr <- which(raw[[1]]=="ProbeSet")[1]; sc <- 4:ncol(raw)
subj<-as.character(raw[3,sc]); visit<-tolower(as.character(raw[2,sc]))
k<-!is.na(subj)&subj!=""; sc<-sc[k]; subj<-subj[k]; visit<-visit[k]
sym<-raw[(hdr+1):nrow(raw),3]
expr<-matrix(as.numeric(as.matrix(raw[(hdr+1):nrow(raw),sc])),nrow=length(sym))
rownames(expr)<-sym; colnames(expr)<-subj; expr<-expr[!is.na(sym)&sym!=""&sym!="---",]
v<-apply(expr,1,var,na.rm=TRUE); o<-order(rownames(expr),-v); expr<-expr[o,]; expr<-expr[!duplicated(rownames(expr)),]

RID<-as.integer(sub(".*_S_","",subj))
dxat<-function(r,vc){x<-am[am$RID==r&am$VISCODE==vc,];if(nrow(x))x$DX[1] else NA}
dxbl<-function(r){x<-am[am$RID==r&am$VISCODE=="bl",];if(nrow(x))x$DX[1] else NA}
g1<-function(r,c){x<-am[am$RID==r,];x<-x[!is.na(x[[c]]),];if(nrow(x))x[[c]][1] else NA}
DX<-mapply(function(r,v){d<-dxat(r,v);if(is.na(d))d<-dxbl(r);d},RID,visit)
grp<-ifelse(DX%in%c("Dementia","AD"),"AD",ifelse(DX=="CN","Control",ifelse(DX=="MCI","MCI",NA)))
AGE<-sapply(RID,g1,"AGE"); SEX<-sapply(RID,function(r)as.character(g1(r,"PTGENDER"))); APOE<-sapply(RID,function(r)suppressWarnings(as.numeric(g1(r,"APOE4"))))

keep<-!is.na(grp)&grp%in%c("Control","AD")
X<-normalizeBetweenArrays(expr[,keep],method="quantile")
s<-read.csv("results/biomarker_survival_after_adjustment.csv",stringsAsFactors=FALSE)
panel<-s$gene[s$tier_orig=="Tier1_Robust"&s$survives%in%c("TRUE","True")]
panel<-intersect(panel,rownames(X))
out<-as.data.frame(t(X[panel,])); out$group<-grp[keep]
out$AGE<-AGE[keep]; out$SEX<-SEX[keep]; out$APOE4<-APOE[keep]
write.csv(out,"results/ai/adni_panel.csv",row.names=FALSE)
cat(sprintf("Saved results/ai/adni_panel.csv: %d samples x %d panel genes (AD=%d, CN=%d)\n",
            nrow(out),length(panel),sum(out$group=="AD"),sum(out$group=="Control")))
