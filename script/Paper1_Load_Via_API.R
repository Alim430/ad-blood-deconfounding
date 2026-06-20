###############################################################################
# Paper 1 — Load Data via GEO API (HTTP, not FTP)
#
# FTP downloads keep timing out. getGEO() uses HTTP which is more reliable.
# Run section by section in RStudio.
###############################################################################

setwd(".")
options(timeout = 3600)  # 1 hour timeout

library(GEOquery)
library(data.table)

raw_dir <- "./data/raw"
dir.create("results", showWarnings = FALSE)

###############################################################################
# SECTION 1: GSE140829 — Parse metadata from truncated series matrix
#
# The expression data (47231 x 2360) is already loaded from raw_data.txt.gz.
# The series matrix is truncated but the METADATA HEADER is intact (first
# ~2749 lines contain all the !Sample_ lines we need).
###############################################################################

cat("=== SECTION 1: GSE140829 Metadata ===\n\n")

# Read ALL lines from the truncated file
con <- gzfile(file.path(raw_dir, "GSE140829_series_matrix.txt.gz"), "r")
lines <- readLines(con, warn = FALSE)
close(con)
cat("Total lines in truncated file:", length(lines), "\n")

# Extract key metadata lines
geo_line     <- grep("^!Sample_geo_accession", lines, value = TRUE)
title_line   <- grep("^!Sample_title", lines, value = TRUE)
diag_lines   <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)

# Parse tab-separated values, removing quotes
parse_meta_line <- function(line) {
  parts <- strsplit(line, "\t")[[1]]
  # Remove the field name (first element)
  parts <- parts[-1]
  # Remove surrounding quotes
  gsub('^"|"$', '', parts)
}

# Get sample IDs
gsm_ids <- parse_meta_line(geo_line)
cat("Number of GSM IDs:", length(gsm_ids), "\n")
cat("First 5:", paste(head(gsm_ids, 5), collapse = ", "), "\n")

# Get titles (contain BeadChip IDs we need for matching)
titles <- parse_meta_line(title_line)
cat("First 3 titles:", paste(head(titles, 3), collapse = "\n  "), "\n")

# Parse characteristics (multiple lines, each is a different characteristic)
cat("\nCharacteristics lines found:", length(diag_lines), "\n")
for (i in seq_along(diag_lines)) {
  vals <- parse_meta_line(diag_lines[i])
  # Show the field name from first value
  cat("  Line", i, ":", substr(vals[1], 1, 60), "  (", length(vals), "values)\n")
}

# Build metadata data frame
meta_140829 <- data.frame(
  geo_accession = gsm_ids,
  title = titles,
  stringsAsFactors = FALSE
)

# Add each characteristics line
for (i in seq_along(diag_lines)) {
  vals <- parse_meta_line(diag_lines[i])
  # Extract field name (before the colon)
  field_name <- sub(":.*", "", vals[1])
  field_name <- trimws(field_name)
  # Extract values (after the colon)
  field_vals <- sub("^[^:]+:\\s*", "", vals)
  meta_140829[[field_name]] <- field_vals
}

cat("\nMeta columns:", paste(colnames(meta_140829), collapse = ", "), "\n")
cat("Dimensions:", nrow(meta_140829), "x", ncol(meta_140829), "\n")

# Show diagnosis distribution
if ("diagnosis" %in% colnames(meta_140829)) {
  cat("\nDiagnosis distribution:\n")
  print(table(meta_140829$diagnosis))
}

# Extract BeadChip ID from title to match with expression data
# Title format: "Whole blood, Control, 200363680054_G [ad_mci]"
meta_140829$beadchip_id <- sub(".*,\\s*([0-9]+_[A-Z])\\s*\\[.*", "\\1", meta_140829$title)
cat("\nFirst 5 BeadChip IDs:", paste(head(meta_140829$beadchip_id, 5), collapse = ", "), "\n")

# Save metadata
write.csv(meta_140829, "results/meta_GSE140829.csv", row.names = FALSE)
cat("Saved: results/meta_GSE140829.csv\n")


###############################################################################
# SECTION 2: GSE63060 — Download via getGEO API (HTTP)
#
# This bypasses FTP entirely and downloads via GEO's web API.
# Should be more reliable than FTP for your connection.
###############################################################################

cat("\n\n=== SECTION 2: GSE63060 via API ===\n")
cat("Downloading via HTTP (not FTP). This may take 10-15 min for ~60 MB...\n")
cat("If it fails, try the browser download link at the end.\n\n")

tryCatch({
  # Force fresh download by removing cached file
  cached <- file.path(raw_dir, "GSE63060_series_matrix.txt.gz")
  if (file.exists(cached) && file.info(cached)$size < 50e6) {
    file.rename(cached, file.path(raw_dir, "GSE63060_series_matrix_truncated.txt.gz"))
    cat("Moved truncated file aside.\n")
  }

  gse63060 <- getGEO("GSE63060", GSEMatrix = TRUE, destdir = raw_dir)
  gse63060 <- gse63060[[1]]
  expr_63060 <- exprs(gse63060)
  meta_63060 <- pData(gse63060)

  cat("SUCCESS!\n")
  cat("Expression:", nrow(expr_63060), "probes x", ncol(expr_63060), "samples\n")

  ch1_cols <- grep(":ch1$", colnames(meta_63060), value = TRUE)
  for (col in ch1_cols) {
    cat("  ", col, ":", paste(head(unique(meta_63060[[col]]), 5), collapse = " | "), "\n")
  }

  write.csv(meta_63060[, c("geo_accession", ch1_cols)],
            "results/meta_GSE63060.csv", row.names = FALSE)
  cat("Saved metadata: results/meta_GSE63060.csv\n")

}, error = function(e) {
  cat("FAILED:", e$message, "\n\n")
  cat("Alternative: Download in Safari/Chrome:\n")
  cat("  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/GSE63060_series_matrix.txt.gz\n")
  cat("  Save to: data/raw/GSE63060_series_matrix.txt.gz\n")
})


###############################################################################
# SECTION 3: GSE63061 — Download via getGEO API
###############################################################################

cat("\n\n=== SECTION 3: GSE63061 via API ===\n")
cat("Downloading via HTTP...\n\n")

tryCatch({
  cached <- file.path(raw_dir, "GSE63061_series_matrix.txt.gz")
  if (file.exists(cached) && file.info(cached)$size < 50e6) {
    file.rename(cached, file.path(raw_dir, "GSE63061_series_matrix_truncated.txt.gz"))
    cat("Moved truncated file aside.\n")
  }

  gse63061 <- getGEO("GSE63061", GSEMatrix = TRUE, destdir = raw_dir)
  gse63061 <- gse63061[[1]]
  expr_63061 <- exprs(gse63061)
  meta_63061 <- pData(gse63061)

  cat("SUCCESS!\n")
  cat("Expression:", nrow(expr_63061), "probes x", ncol(expr_63061), "samples\n")

  ch1_cols <- grep(":ch1$", colnames(meta_63061), value = TRUE)
  for (col in ch1_cols) {
    cat("  ", col, ":", paste(head(unique(meta_63061[[col]]), 5), collapse = " | "), "\n")
  }

  write.csv(meta_63061[, c("geo_accession", ch1_cols)],
            "results/meta_GSE63061.csv", row.names = FALSE)

}, error = function(e) {
  cat("FAILED:", e$message, "\n\n")
  cat("Alternative: Download in Safari/Chrome:\n")
  cat("  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/GSE63061_series_matrix.txt.gz\n")
})


###############################################################################
# SECTION 4: GSE85426 — Get expression data
#
# The series matrix has metadata (180 samples) but 0 expression probes.
# Expression might be in supplementary files or need non-normalized data.
###############################################################################

cat("\n\n=== SECTION 4: GSE85426 Expression Data ===\n")

# Check if we already have metadata from previous run
if (!exists("meta_85426")) {
  cat("Loading metadata...\n")
  gse85426 <- getGEO("GSE85426", GSEMatrix = TRUE, destdir = raw_dir)
  gse85426 <- gse85426[[1]]
  meta_85426 <- pData(gse85426)
}
cat("Metadata:", nrow(meta_85426), "samples\n")

# Check supplementary files for this dataset
cat("\nChecking for supplementary expression files...\n")
cat("This dataset (Tan et al., Malaysian cohort) uses Agilent microarray.\n")
cat("The expression data may be in the supplementary files on GEO.\n\n")

# Try to get supplementary file info
tryCatch({
  gse85426_full <- getGEO("GSE85426", GSEMatrix = FALSE)
  # Check what supplementary files exist
  supp <- Meta(gse85426_full)$supplementary_file
  cat("Supplementary files listed in GEO record:\n")
  if (length(supp) > 0) {
    for (s in supp) cat("  ", s, "\n")
  } else {
    cat("  (none listed at study level)\n")
  }

  # Check first sample for supplementary files
  first_gsm <- GSMList(gse85426_full)[[1]]
  gsm_supp <- Meta(first_gsm)$supplementary_file
  cat("\nFirst sample supplementary files:\n")
  if (length(gsm_supp) > 0) {
    for (s in gsm_supp) cat("  ", s, "\n")
  } else {
    cat("  (none)\n")
  }
}, error = function(e) {
  cat("Error checking supplementary files:", e$message, "\n")
})

# Save metadata regardless
ch1_cols <- grep(":ch1$", colnames(meta_85426), value = TRUE)
write.csv(meta_85426[, c("geo_accession", ch1_cols)],
          "results/meta_GSE85426.csv", row.names = FALSE)
cat("\nSaved metadata: results/meta_GSE85426.csv\n")

cat("\nNOTE: If GSE85426 has no accessible expression data, we can still\n")
cat("proceed with the other 4 datasets (~2800+ samples). GSE85426 (180 samples)\n")
cat("adds the Malaysian/Asian cohort for cross-ethnicity analysis, but is not\n")
cat("strictly required for the core RS/FM pipeline.\n")


###############################################################################
# SECTION 5: GSE270454 — Prepare metadata from column names
#
# Column names already encode condition: MCI_1, AD_1, ASO_1, ASM_1
# ASO = Asymptomatic Obese, ASM = Asymptomatic Metabolic (check paper)
###############################################################################

cat("\n\n=== SECTION 5: GSE270454 Metadata from Column Names ===\n")

f270 <- file.path(raw_dir, "GSE270454_RNAseq-combined-counts-matrix.csv.gz")
counts_270454 <- fread(f270, data.table = FALSE)
cat("Full count matrix:", nrow(counts_270454), "genes x", ncol(counts_270454)-1, "samples\n")

# Gene names are in V1
gene_names <- counts_270454$V1
counts_270454$V1 <- NULL
sample_names <- colnames(counts_270454)

# Extract condition from column names
conditions <- sub("_[0-9]+$", "", sample_names)
cat("\nCondition distribution:\n")
print(table(conditions))

meta_270454 <- data.frame(
  sample_id = sample_names,
  condition = conditions,
  stringsAsFactors = FALSE
)

# Also get metadata from series matrix if available
tryCatch({
  gse270454 <- getGEO("GSE270454", GSEMatrix = TRUE, destdir = raw_dir)
  gse270454_data <- gse270454[[1]]
  meta_270454_full <- pData(gse270454_data)
  cat("\nFull metadata from GEO:", nrow(meta_270454_full), "samples\n")
  ch1_cols <- grep(":ch1$", colnames(meta_270454_full), value = TRUE)
  for (col in ch1_cols) {
    cat("  ", col, ":", paste(head(unique(meta_270454_full[[col]]), 6), collapse = " | "), "\n")
  }
  write.csv(meta_270454_full[, c("geo_accession", "title", ch1_cols)],
            "results/meta_GSE270454.csv", row.names = FALSE)
}, error = function(e) {
  cat("Could not get full metadata:", e$message, "\n")
  write.csv(meta_270454, "results/meta_GSE270454.csv", row.names = FALSE)
})


###############################################################################
# SECTION 6: GSE248423 — Check what's available
###############################################################################

cat("\n\n=== SECTION 6: GSE248423 ===\n")

tryCatch({
  gse248423 <- getGEO("GSE248423", GSEMatrix = TRUE, destdir = raw_dir)
  gse248423_data <- gse248423[[1]]
  expr_248423 <- exprs(gse248423_data)
  meta_248423 <- pData(gse248423_data)

  cat("Expression:", nrow(expr_248423), "x", ncol(expr_248423), "\n")
  cat("Metadata:", nrow(meta_248423), "samples\n")

  ch1_cols <- grep(":ch1$", colnames(meta_248423), value = TRUE)
  for (col in ch1_cols) {
    cat("  ", col, ":", paste(head(unique(meta_248423[[col]]), 6), collapse = " | "), "\n")
  }

  # Check supplementary files
  cat("\nChecking supplementary files...\n")
  gse248423_raw <- getGEO("GSE248423", GSEMatrix = FALSE)
  supp <- Meta(gse248423_raw)$supplementary_file
  if (length(supp) > 0) {
    cat("Supplementary files:\n")
    for (s in supp) cat("  ", s, "\n")
  }

  write.csv(meta_248423[, c("geo_accession", ch1_cols)],
            "results/meta_GSE248423.csv", row.names = FALSE)
}, error = function(e) {
  cat("Error:", e$message, "\n")
})


###############################################################################
# SECTION 7: FINAL STATUS
###############################################################################

cat("\n\n")
cat("================================================================\n")
cat("  FINAL DATA LOADING STATUS\n")
cat("================================================================\n\n")

status <- data.frame(
  Dataset = c("GSE140829", "GSE63060", "GSE63061", "GSE85426", "GSE270454", "GSE248423"),
  Type = c("Microarray", "Microarray", "Microarray", "Microarray", "RNA-seq", "RNA-seq"),
  Expression = c(
    ifelse(exists("counts_140829"), "LOADED", "MISSING"),
    ifelse(exists("expr_63060") && nrow(expr_63060) > 0, "LOADED", "MISSING"),
    ifelse(exists("expr_63061") && nrow(expr_63061) > 0, "LOADED", "MISSING"),
    ifelse(exists("expr_85426") && nrow(expr_85426) > 0, "LOADED", "MISSING"),
    ifelse(exists("counts_270454"), "LOADED", "MISSING"),
    ifelse(exists("expr_248423") && nrow(expr_248423) > 0, "LOADED", "MISSING")
  ),
  Metadata = c(
    ifelse(exists("meta_140829"), "LOADED", "MISSING"),
    ifelse(exists("meta_63060"), "LOADED", "MISSING"),
    ifelse(exists("meta_63061"), "LOADED", "MISSING"),
    ifelse(exists("meta_85426"), "LOADED", "MISSING"),
    "FROM_COLNAMES",
    ifelse(exists("meta_248423"), "LOADED", "MISSING")
  ),
  stringsAsFactors = FALSE
)

print(status)

cat("\nSaving workspace...\n")
save.image("results/01_all_loaded_data.RData")
cat("Done! Saved to results/01_all_loaded_data.RData\n")
cat("\nNext step: Run Paper1_R_Pipeline.R starting from Section 4 (preprocessing)\n")
