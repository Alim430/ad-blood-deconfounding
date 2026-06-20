###############################################################################
# AI Module 2a: export scRNA counts for deep embedding  [SERVER, GPU MODE]
# Writes a 10x-style mtx so Python (scVI / scGPT) can read it without SeuratDisk.
#   Rscript script/ai/Step20a_export_h5.R
###############################################################################
setwd(Sys.getenv("AD_PAPER_ROOT", "."))
suppressWarnings(suppressMessages({library(Seurat); library(Matrix)}))
dir.create("results/ai", recursive=TRUE, showWarnings=FALSE)
# auto-find the processed Seurat object (name may vary)
rds <- Sys.getenv("RDS","")
if (rds=="") { cand <- list.files("results", pattern="\\.rds$", full.names=TRUE)
  pref <- cand[grepl("seurat|processed|GSE2266|clustered", cand, ignore.case=TRUE)]
  rds <- if (length(pref)) pref[which.max(file.size(pref))] else cand[which.max(file.size(cand))] }
cat("Using Seurat object:", rds, "\n")
pbmc <- readRDS(rds)
# export only highly-variable genes (scVI standard) — keeps the mtx ~10x smaller
vf <- VariableFeatures(pbmc)
if (length(vf) < 100) { pbmc <- FindVariableFeatures(pbmc, nfeatures=2000); vf <- VariableFeatures(pbmc) }
cat("Exporting", length(vf), "highly-variable genes\n")
cnt  <- GetAssayData(pbmc, assay="RNA", layer="counts")[vf, , drop=FALSE]
writeMM(cnt, "results/ai/counts.mtx")
writeLines(rownames(cnt), "results/ai/genes.txt")
write.csv(data.frame(barcode=colnames(cnt), condition=pbmc$condition,
                     cell_type=pbmc$cell_type),
          "results/ai/cells.csv", row.names=FALSE)
cat("Saved results/ai/{counts.mtx, genes.txt, cells.csv} — feed to Step20b.\n")
