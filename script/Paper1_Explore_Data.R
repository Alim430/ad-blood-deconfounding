###############################################################################
# Paper 1 — Explore What We Have So Far
#
# GSE140829 loaded successfully (47231 probes x 2360 samples)
# BUT: it's Illumina BeadChip MICROARRAY, not RNA-seq!
#   Evidence: ILMN_* probe IDs, columns = PROBE_ID, SYMBOL, etc.
#   This means: process with limma, not DESeq2
#
# Run this to understand the data structure before the series matrix
# re-downloads finish.
###############################################################################

library(data.table)

setwd(".")
raw_dir <- "./data/raw"

# ========================================================================
# PART 1: Explore GSE140829 structure
# ========================================================================

cat("=== GSE140829 RAW DATA STRUCTURE ===\n\n")

# Reload if not in memory
if (!exists("counts_140829_raw")) {
  cat("Loading GSE140829_raw_data.txt.gz...\n")
  counts_140829_raw <- fread(file.path(raw_dir, "GSE140829_raw_data.txt.gz"),
                             data.table = FALSE, nrows = 10)
}

# Show ALL column names
cat("Total columns:", ncol(counts_140829_raw), "\n\n")
cat("First 20 column names:\n")
print(colnames(counts_140829_raw)[1:20])

cat("\nLast 10 column names:\n")
print(colnames(counts_140829_raw)[(ncol(counts_140829_raw)-9):ncol(counts_140829_raw)])

# The annotation columns (non-sample) are things like:
# PROBE_ID, SYMBOL, SEARCH_KEY, ILMN_GENE, CHROMOSOME, etc.
# Sample columns will be GSM IDs or patient IDs

# Identify which columns are annotation vs expression
cat("\n--- Identifying annotation vs sample columns ---\n")
# Annotation columns usually have text, sample columns have numbers
first_row <- counts_140829_raw[1, ]
annotation_cols <- c()
sample_cols <- c()
for (col in colnames(counts_140829_raw)) {
  val <- first_row[[col]]
  if (is.numeric(val) || suppressWarnings(!is.na(as.numeric(val)))) {
    sample_cols <- c(sample_cols, col)
  } else {
    annotation_cols <- c(annotation_cols, col)
  }
}

cat("Annotation columns:", length(annotation_cols), "\n")
print(annotation_cols)
cat("\nSample columns (first 10 of", length(sample_cols), "):\n")
print(head(sample_cols, 10))

# Check if sample columns look like GSM IDs
cat("\nDo sample columns look like GSM IDs?\n")
gsm_pattern <- grep("^GSM", sample_cols, value = TRUE)
cat("  GSM-pattern matches:", length(gsm_pattern), "\n")
if (length(gsm_pattern) > 0) {
  cat("  Examples:", paste(head(gsm_pattern, 5), collapse = ", "), "\n")
} else {
  cat("  Sample naming pattern:", paste(head(sample_cols, 5), collapse = ", "), "\n")
}


# ========================================================================
# PART 2: Try to get metadata from truncated series matrix
# ========================================================================

cat("\n\n=== TRYING TO EXTRACT METADATA FROM TRUNCATED FILE ===\n")

meta_file <- file.path(raw_dir, "GSE140829_series_matrix.txt.gz")
if (file.exists(meta_file)) {
  cat("Reading lines from truncated series matrix...\n")
  con <- gzfile(meta_file, "r")
  lines <- readLines(con, n = 5000)  # Read first 5000 lines
  close(con)
  cat("Read", length(lines), "lines\n")

  # Extract sample characteristics
  char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  title_lines <- grep("^!Sample_title", lines, value = TRUE)
  geo_lines <- grep("^!Sample_geo_accession", lines, value = TRUE)
  source_lines <- grep("^!Sample_source_name_ch1", lines, value = TRUE)

  cat("\nMetadata lines found:\n")
  cat("  Sample titles:", length(title_lines), "\n")
  cat("  GEO accessions:", length(geo_lines), "\n")
  cat("  Characteristics:", length(char_lines), "\n")
  cat("  Source names:", length(source_lines), "\n")

  # Show what characteristics look like
  if (length(char_lines) > 0) {
    cat("\nFirst characteristics line (truncated):\n")
    cat(substr(char_lines[1], 1, 200), "\n...\n")
  }

  if (length(title_lines) > 0) {
    cat("\nFirst title line (truncated):\n")
    cat(substr(title_lines[1], 1, 200), "\n...\n")
  }

  if (length(source_lines) > 0) {
    cat("\nSource name line (truncated):\n")
    cat(substr(source_lines[1], 1, 200), "\n...\n")
  }
}


# ========================================================================
# PART 3: Explore GSE270454 count matrix
# ========================================================================

cat("\n\n=== GSE270454 COUNT MATRIX ===\n")
f270 <- file.path(raw_dir, "GSE270454_RNAseq-combined-counts-matrix.csv.gz")
if (file.exists(f270)) {
  counts_270454 <- fread(f270, data.table = FALSE, nrows = 5)
  cat("Dimensions (first 5 rows):", ncol(counts_270454), "columns\n")
  cat("Column names:\n")
  print(colnames(counts_270454))
  cat("\nFirst 5 gene IDs:", paste(counts_270454[[1]][1:5], collapse = ", "), "\n")
}


# ========================================================================
# PART 4: Try loading GSE85426 via GEO API (small dataset)
# ========================================================================

cat("\n\n=== GSE85426 — Trying GEO API download ===\n")
cat("This downloads the full dataset directly from GEO (~180 samples)...\n")

library(GEOquery)
tryCatch({
  gse85426 <- getGEO("GSE85426", GSEMatrix = TRUE, destdir = raw_dir)
  gse85426 <- gse85426[[1]]
  expr_85426 <- exprs(gse85426)
  meta_85426 <- pData(gse85426)
  cat("SUCCESS!\n")
  cat("Expression:", nrow(expr_85426), "probes x", ncol(expr_85426), "samples\n")
  cat("Metadata:", nrow(meta_85426), "samples\n")

  ch1_cols <- grep(":ch1$", colnames(meta_85426), value = TRUE)
  cat("\nMetadata columns:\n")
  for (col in ch1_cols) {
    vals <- unique(meta_85426[[col]])
    cat("  ", col, ":", paste(head(vals, 6), collapse = " | "), "\n")
  }
}, error = function(e) {
  cat("Failed:", e$message, "\n")
  cat("The series matrix for GSE85426 might genuinely be metadata-only.\n")
  cat("Try downloading via getGEO with GSEMatrix=FALSE and check.\n")
})


# ========================================================================
# PART 5: Summary
# ========================================================================

cat("\n\n=== DATASET STATUS SUMMARY ===\n\n")
cat("IMPORTANT CORRECTION: GSE140829 is MICROARRAY (Illumina BeadChip),\n")
cat("not RNA-seq! The ILMN_* probe IDs confirm this.\n\n")

cat("Dataset Classification:\n")
cat("  MICROARRAY (process with limma):\n")
cat("    - GSE140829: Illumina BeadChip, 47231 probes x 2360 samples [LOADED]\n")
cat("    - GSE63060:  Illumina HumanHT-12, ~329 samples [NEED SERIES MATRIX]\n")
cat("    - GSE63061:  Illumina HumanHT-12, ~382 samples [NEED SERIES MATRIX]\n")
cat("    - GSE85426:  Agilent, ~180 samples [TRYING API]\n")
cat("  RNA-SEQ (process with DESeq2):\n")
cat("    - GSE270454: RNA-seq, ~45 samples [COUNT MATRIX LOADED]\n")
cat("    - GSE248423: RNA-seq, ~100 samples [NEED COUNT FILE]\n\n")

cat("Blocked on:\n")
cat("  1. Re-download GSE140829 series matrix (metadata) — 216 MB\n")
cat("  2. Re-download GSE63060 series matrix — 60 MB\n")
cat("  3. Re-download GSE63061 series matrix — 60 MB\n")
cat("  4. Find GSE248423 count file\n")
cat("\nRun fix_downloads_v2.sh in Terminal to fix #1-3.\n")
cat("For #4, visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248423\n")
