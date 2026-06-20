###############################################################################
# Step 22b [SERVER]: aggregate FULL-GENE pseudobulk (donor x cell-type) ONLY.
# No DE here (scrna env lacks DESeq2/edgeR). Saves a tiny pseudobulk object;
# run the DE on the Mac with Step22c. Needs the 90GB instance (loads 7GB object).
#   Rscript script/ai/Step22b_server_pseudobulk.R
# Output: results/pseudobulk_fullgene.rds   (download to Mac)
###############################################################################
setwd(Sys.getenv("AD_PAPER_ROOT", "."))
suppressWarnings(suppressMessages({library(Seurat); library(Matrix)}))
rds <- "results/GSE226602_seurat_processed.rds"
cat("Loading", rds, "(one time)...\n"); obj <- readRDS(rds)
cnt <- GetAssayData(obj, assay="RNA", layer="counts")
md  <- data.frame(barcode=colnames(obj), cell_type=as.character(obj$cell_type),
                  condition=as.character(obj$condition), stringsAsFactors=FALSE)
md$donor <- sub("_[ACGT]+-1$","",md$barcode)
md <- md[md$condition %in% c("AD","Control") & !is.na(md$cell_type), ]
cnt <- cnt[, md$barcode, drop=FALSE]
md$group <- paste(md$donor, md$cell_type, sep="__")
cat(sprintf("%d genes x %d cells; %d donors\n", nrow(cnt), ncol(cnt), length(unique(md$donor))))

grp <- factor(md$group)
ind <- sparse.model.matrix(~0+grp); colnames(ind) <- levels(grp)
pb  <- as.matrix(cnt %*% ind); storage.mode(pb) <- "integer"     # genes x groups
gm  <- unique(md[,c("group","donor","cell_type","condition")]); rownames(gm) <- gm$group
gm  <- gm[colnames(pb),]; gm$ncell <- as.integer(table(md$group)[colnames(pb)])
rownames(pb) <- rownames(cnt)
saveRDS(list(pb=pb, gm=gm), "results/pseudobulk_fullgene.rds")
cat(sprintf("Saved results/pseudobulk_fullgene.rds  (%d genes x %d pseudobulk samples)\n",
            nrow(pb), ncol(pb)))
cat("Download to Mac, then run: Rscript script/ai/Step22c_mac_pseudobulk_DE.R\n")
