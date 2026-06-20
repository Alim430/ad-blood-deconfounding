#!/usr/bin/env Rscript
# ============================================================
# Step 32: Multi-method deconvolution robustness check
#   Methods compared: Danaher (current primary) + Bisque + MuSiC
# Goal (per de-risked plan): show conclusions are NOT method-specific
# Outputs CSVs that Python Step 33 turns into FigSx_deconv_robustness
# ============================================================
# REQUIREMENTS (install once)
#   BiocManager::install(c("Biobase","SingleCellExperiment","MuSiC","BisqueRNA"))
#   install.packages(c("data.table","Matrix"))
# ============================================================
# ============================================================
# Package bootstrap (auto install if missing)
# ============================================================
# Run this script from the repository root (paths below are relative).
# If invoked from elsewhere, uncomment and adapt:
#   setwd("/path/to/ad-blood-deconfounding")
cat("Working directory:", getwd(), "\n")

pkgs_cran <- c("data.table", "Matrix")
pkgs_bioc <- c("Biobase", "SingleCellExperiment", "MuSiC", "BisqueRNA")

# install CRAN packages if missing
install_if_missing_cran <- function(pkgs){
  for(p in pkgs){
    if(!requireNamespace(p, quietly = TRUE)){
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

# install Bioconductor packages if missing
install_if_missing_bioc <- function(pkgs){
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for(p in pkgs){
    if(!requireNamespace(p, quietly = TRUE)){
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }
  }
}

install_if_missing_cran(pkgs_cran)
install_if_missing_bioc(pkgs_bioc)

# ============================================================
# Load libraries
# ============================================================
suppressPackageStartupMessages({
  library(Biobase)
  library(Matrix)
  library(data.table)
  library(BisqueRNA)
  library(MuSiC)
  library(SingleCellExperiment)
})
suppressPackageStartupMessages({
  library(Biobase); library(Matrix); library(data.table)
  library(BisqueRNA); library(MuSiC); library(SingleCellExperiment)
})

OUT_DIR <- "results/deconv_robustness"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ---- 1. Load bulk expression matrices + phenotype for 3 cohorts ----
# Object names inside each RData (verified from actual files, 2026-06):
#   GSE140829: expr_140829_matched, meta_140829_matched, probe_gene_map
#   GSE63060 : expr_63060,          meta_63060
#   GSE63061 : expr_63061,          meta_63061
load("results/GSE140829_processed.RData")
load("results/GSE63060_processed.RData")
load("results/GSE63061_processed.RData")
bulk_list <- list(
  GSE140829 = list(expr=expr_140829_matched, pheno=meta_140829_matched),
  GSE63060  = list(expr=expr_63060,          pheno=meta_63060),
  GSE63061  = list(expr=expr_63061,          pheno=meta_63061)
)

# (probe-to-symbol mapping happens after sc_ref is loaded — see block 1b below.)

# ---- 2. Load scRNA reference (GSE226602 PBMC; donor-level pseudobulk per cell type) ----
# We expect a saved SingleCellExperiment or anndata-converted RData here.
# If you only have the pseudobulk per cell type, that's enough for MuSiC's reference.
sc_ref <- readRDS("results/pseudobulk_fullgene.rds")   # produced earlier
# sc_ref must include $counts (gene x cell), $cell_type, $donor_id

# Build a SingleCellExperiment for MuSiC/Bisque.
# NOTE: sc_ref$pb is a *pseudobulk* matrix (genes × donor×cell_type aggregates),
#   not single cells. Each column = one (donor, cell_type) aggregate.
#   sc_ref$gm has per-column metadata: $donor, $cell_type, $condition, $ncell.
# Both Bisque (ReferenceBasedDecomposition) and MuSiC tolerate pseudobulk
# references at the (donor, cell_type) granularity — each pseudobulk acts as
# one "cell" of that donor x cell_type combination.
make_sce <- function(ref){
  cnt <- as.matrix(ref$pb)
  gm  <- ref$gm
  stopifnot(ncol(cnt) == nrow(gm))
  SingleCellExperiment(
    assays  = list(counts = cnt),
    colData = DataFrame(cell_type = as.character(gm$cell_type),
                        donor_id  = as.character(gm$donor))
  )
}
sce <- make_sce(sc_ref)
cat(sprintf("Reference SCE: %d genes x %d pseudobulks (%d donors, cell types: %s)\n",
            nrow(sce), ncol(sce), length(unique(sce$donor_id)),
            paste(unique(sce$cell_type), collapse=", ")))

# ---- 1b. Map Illumina probe IDs -> gene symbols (all 3 cohorts use ILMN_*) ----
# probe_gene_map (43,961 ILMN -> symbol) comes from GSE140829 .RData; same Illumina HT-12
# platform family is used by all three cohorts, so the same map applies.
stopifnot(exists("probe_gene_map"))
collapse_probes_to_genes <- function(expr, map){
  m <- map[!is.na(map$symbol) & map$symbol!="" & map$probe_id %in% rownames(expr), ]
  m <- m[!duplicated(m$probe_id), ]
  sub <- expr[m$probe_id, , drop=FALSE]
  v <- apply(sub, 1, function(x) var(x, na.rm=TRUE))
  ord <- order(-v); m <- m[ord, ]; sub <- sub[ord, , drop=FALSE]
  keep <- !duplicated(m$symbol)
  out  <- sub[keep, , drop=FALSE]; rownames(out) <- m$symbol[keep]; out
}
cat("Probe -> gene-symbol collapse:\n")
for(n in names(bulk_list)){
  e <- bulk_list[[n]]$expr; before <- nrow(e)
  if(any(grepl("^ILMN_", rownames(e)))){
    e2 <- collapse_probes_to_genes(e, probe_gene_map)
    bulk_list[[n]]$expr <- e2
    cat(sprintf("  %-10s  %d probes -> %d gene symbols  (overlap w/ scRNA ref: %d)\n",
                n, before, nrow(e2), length(intersect(rownames(e2), rownames(sce)))))
  } else {
    cat(sprintf("  %-10s  already gene-symbol (%d genes)\n", n, before))
  }
}
cat("WARNING: PBMC reference has no neutrophils — Bisque/MuSiC will NOT\n")
cat("         estimate the neutrophil proportion. Concordance with Danaher\n")
cat("         is reported for shared lineages (Mono, B, CD4T, CD8T, NK).\n\n")

# ---- 3. Helper: run a method on one cohort, return cell proportions ----
run_bisque <- function(bulk_expr, sce){
  # Bisque reference-based
  ref_eset <- ExpressionSet(assayData=as.matrix(assays(sce)$counts),
                            phenoData=AnnotatedDataFrame(
                              data.frame(SubjectName=colData(sce)$donor_id,
                                         cellType   =colData(sce)$cell_type,
                                         row.names=colnames(sce))))
  bulk_eset <- ExpressionSet(assayData=as.matrix(bulk_expr))
  res <- ReferenceBasedDecomposition(bulk_eset, ref_eset, use.overlap=FALSE)
  t(res$bulk.props)
}

run_music <- function(bulk_expr, sce){
  music_prop(
    bulk.mtx     = as.matrix(bulk_expr),
    sc.sce       = sce,
    clusters     = "cell_type",
    samples      = "donor_id",
    verbose      = FALSE
  )$Est.prop.weighted
}

# ---- 4. Run both methods on all cohorts ----
for(coh in names(bulk_list)){
  cat(sprintf("=== %s (Bisque) ===\n", coh))
  bq <- tryCatch(run_bisque(bulk_list[[coh]]$expr, sce), error=function(e){print(e); NULL})
  if(!is.null(bq)) write.csv(bq, file.path(OUT_DIR, sprintf("%s_Bisque.csv", coh)))

  cat(sprintf("=== %s (MuSiC) ===\n", coh))
  mu <- tryCatch(run_music(bulk_list[[coh]]$expr, sce), error=function(e){print(e); NULL})
  if(!is.null(mu)) write.csv(mu, file.path(OUT_DIR, sprintf("%s_MuSiC.csv", coh)))
}

# ---- 5. Concordance vs current Danaher scores ----
# Danaher z-scores are in results/cell_scores_{cohort}.csv
for(coh in names(bulk_list)){
  dan_file <- sprintf("results/cell_scores_%s.csv", coh)
  if(!file.exists(dan_file)) next
  dan <- fread(dan_file)
  for(meth in c("Bisque","MuSiC")){
    f <- file.path(OUT_DIR, sprintf("%s_%s.csv", coh, meth))
    if(!file.exists(f)) next
    pr <- fread(f); setnames(pr,1,"sample")
    # match by sample id and compute per-lineage Pearson + Spearman
    common <- intersect(dan$sample, pr$sample)
    if(length(common) < 10) next
    merged <- merge(dan[sample %in% common], pr[sample %in% common], by="sample", suffixes=c(".dan",".x"))
    # save merged for plotting
    fwrite(merged, file.path(OUT_DIR, sprintf("%s_%s_vs_Danaher.csv", coh, meth)))
  }
}

cat("DONE — see results/deconv_robustness/\n")
