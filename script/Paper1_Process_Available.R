###############################################################################
# Paper 1 — Process What We Have + Extract GSE248423 from SOFT
#
# Status:
#   GSE140829: Expression + Metadata LOADED (587 samples)
#   GSE270454: Expression + Metadata LOADED (45 samples)
#   GSE248423: SOFT data downloaded — extract expression from per-sample tables
#   GSE85426:  Metadata only (drop for now, 180 samples — not critical)
#   GSE63060:  BLOCKED — download in Safari (see bottom of script)
#   GSE63061:  BLOCKED — download in Safari (see bottom of script)
#
# Run in RStudio section by section.
###############################################################################

setwd(".")
options(timeout = 3600)

library(GEOquery)
library(data.table)
library(limma)

raw_dir <- "./data/raw"
dir.create("results", showWarnings = FALSE)


###############################################################################
# SECTION 1: Process GSE140829 Expression (extract AVG_SIGNAL only)
#
# The raw data has 4 columns per sample:
#   {BeadChipID}.AVG_SIGNAL    <- this is what we need
#   {BeadChipID}.DETECTION PVAL
#   {BeadChipID}.BEAD_STDERR
#   {BeadChipID}.AVG_NBEADS
#
# Plus 12 annotation columns (PROBE_ID, SYMBOL, etc.)
###############################################################################

cat("=== SECTION 1: GSE140829 — Extract Expression Matrix ===\n\n")

# Load full data (only if not already in memory)
if (!exists("gse140829_raw") || !is.data.frame(gse140829_raw)) {
  cat("Loading full GSE140829 raw data (327 MB, ~2 min)...\n")
  gse140829_raw <- fread(file.path(raw_dir, "GSE140829_raw_data.txt.gz"),
                         data.table = FALSE)
  cat("Loaded:", nrow(gse140829_raw), "probes x", ncol(gse140829_raw), "columns\n")
}

# Identify annotation columns
annotation_cols <- c("PROBE_ID", "SYMBOL", "SEARCH_KEY", "ILMN_GENE",
                     "CHROMOSOME", "DEFINITION", "SYNONYMS", "TRANSCRIPT",
                     "SOURCE_REFERENCE_ID", "REFSEQ_ID", "UNIGENE_ID",
                     "ACCESSION", "SOURCE")

# Extract only AVG_SIGNAL columns (case-insensitive match)
all_cols <- colnames(gse140829_raw)
signal_cols <- grep("AVG_Signal|AVG_SIGNAL", all_cols, value = TRUE, ignore.case = TRUE)
cat("Found", length(signal_cols), "AVG_SIGNAL columns (= samples)\n")

# Also grab detection p-value columns for filtering
pval_cols <- grep("DETECTION.PVAL|Detection.Pval", all_cols, value = TRUE, ignore.case = TRUE)
cat("Found", length(pval_cols), "detection p-value columns\n")

# Extract signal matrix
probe_ids <- gse140829_raw$PROBE_ID
symbols <- gse140829_raw$SYMBOL
expr_140829 <- as.matrix(gse140829_raw[, signal_cols])
rownames(expr_140829) <- probe_ids
cat("Expression matrix:", nrow(expr_140829), "probes x", ncol(expr_140829), "samples\n")

# Clean up column names to get BeadChip IDs
# "5872617031_A.AVG_SIGNAL" -> "5872617031_A"
sample_ids_140829 <- sub("\\.AVG_Signal$|\\.AVG_SIGNAL$", "", colnames(expr_140829),
                         ignore.case = TRUE)
colnames(expr_140829) <- sample_ids_140829
cat("Sample ID examples:", paste(head(sample_ids_140829, 3), collapse = ", "), "\n")

# Load metadata (from Section 1 of previous script)
if (!exists("meta_140829") || !"beadchip_id" %in% colnames(meta_140829)) {
  meta_140829 <- read.csv("results/meta_GSE140829.csv", stringsAsFactors = FALSE)
}
cat("Metadata:", nrow(meta_140829), "samples\n")

# Match expression columns to metadata via BeadChip ID
matched <- meta_140829$beadchip_id %in% sample_ids_140829
cat("Metadata samples matched to expression:", sum(matched), "of", nrow(meta_140829), "\n")

if (sum(matched) < nrow(meta_140829)) {
  # Check what's not matching
  unmatched <- meta_140829$beadchip_id[!matched]
  cat("Unmatched BeadChip IDs (first 5):", paste(head(unmatched, 5), collapse = ", "), "\n")
  cat("Expression column IDs (first 5):", paste(head(sample_ids_140829, 5), collapse = ", "), "\n")
}

# Subset expression to matched samples only, in metadata order
expr_140829_matched <- expr_140829[, meta_140829$beadchip_id[matched], drop = FALSE]
meta_140829_matched <- meta_140829[matched, ]
cat("\nFinal matched: ", ncol(expr_140829_matched), "samples\n")
cat("Diagnosis distribution in matched data:\n")
print(table(meta_140829_matched$diagnosis))

# Probe-to-gene mapping
probe_gene_map <- data.frame(
  probe_id = probe_ids,
  symbol = symbols,
  stringsAsFactors = FALSE
)
# Remove probes with no gene symbol
probe_gene_map <- probe_gene_map[!is.na(probe_gene_map$symbol) &
                                  probe_gene_map$symbol != "", ]
cat("\nProbes with gene symbols:", nrow(probe_gene_map), "of", length(probe_ids), "\n")

# Save processed data
save(expr_140829_matched, meta_140829_matched, probe_gene_map,
     file = "results/GSE140829_processed.RData")
cat("Saved: results/GSE140829_processed.RData\n")


###############################################################################
# SECTION 2: Process GSE270454 (RNA-seq counts)
###############################################################################

cat("\n\n=== SECTION 2: GSE270454 — Process RNA-seq Counts ===\n\n")

f270 <- file.path(raw_dir, "GSE270454_RNAseq-combined-counts-matrix.csv.gz")
counts_270454 <- fread(f270, data.table = FALSE)
gene_names_270454 <- counts_270454$V1
counts_270454$V1 <- NULL
rownames(counts_270454) <- gene_names_270454
counts_270454 <- as.matrix(counts_270454)

cat("Count matrix:", nrow(counts_270454), "genes x", ncol(counts_270454), "samples\n")

# Build metadata from column names
meta_270454 <- data.frame(
  sample_id = colnames(counts_270454),
  condition = sub("_[0-9]+$", "", colnames(counts_270454)),
  stringsAsFactors = FALSE
)

# For AD vs Control analysis, we'll use AD and MCI
# ASO = Asymptomatic Obese, ASM = Asymptomatic Metabolic — need to check paper
# For now, keep all conditions
cat("Condition distribution:\n")
print(table(meta_270454$condition))

# Note: This dataset has NO healthy controls!
# AD=10, ASM=11, ASO=14, MCI=10
# ASM and ASO might serve as "non-demented" comparisons
# Need to verify from the paper what these abbreviations mean
cat("\nWARNING: No 'Control' group in this dataset.\n")
cat("ASO/ASM may be non-demented comparison groups — verify from paper.\n")

save(counts_270454, meta_270454, file = "results/GSE270454_processed.RData")
cat("Saved: results/GSE270454_processed.RData\n")


###############################################################################
# SECTION 3: Extract GSE248423 expression from SOFT data
#
# The getGEO(GSEMatrix=FALSE) call downloaded all 197 GSM entities.
# Each GSM may have a data table with expression values.
###############################################################################

cat("\n\n=== SECTION 3: GSE248423 — Extract from SOFT ===\n\n")

# Check if gse248423_raw object exists from previous run
if (!exists("gse248423_raw")) {
  cat("Re-downloading GSE248423 SOFT data (this worked before)...\n")
  tryCatch({
    gse248423_raw <- getGEO("GSE248423", GSEMatrix = FALSE)
    cat("Downloaded.\n")
  }, error = function(e) {
    cat("Download failed:", e$message, "\n")
    cat("Skipping GSE248423 for now.\n")
  })
}

if (exists("gse248423_raw")) {
  # Check if individual GSMs have data tables
  gsm_list <- GSMList(gse248423_raw)
  cat("Number of GSMs:", length(gsm_list), "\n")

  # Check first GSM for data table
  first_gsm <- gsm_list[[1]]
  dt <- Table(first_gsm)
  cat("First GSM data table dimensions:", nrow(dt), "x", ncol(dt), "\n")

  if (nrow(dt) > 0) {
    cat("Columns:", paste(colnames(dt), collapse = ", "), "\n")
    cat("First 3 rows:\n")
    print(head(dt, 3))

    # Extract expression from all GSMs
    cat("\nExtracting expression from all GSMs...\n")
    all_gsm_names <- names(gsm_list)

    # Determine which column has the expression value
    # Common names: VALUE, SIGNAL, COUNT, TPM, FPKM
    value_col <- intersect(colnames(dt), c("VALUE", "SIGNAL", "COUNT",
                                            "TPM", "FPKM", "COUNTS"))
    if (length(value_col) == 0) {
      # Try numeric columns
      num_cols <- sapply(dt, is.numeric)
      value_col <- colnames(dt)[num_cols]
      cat("Numeric columns found:", paste(value_col, collapse = ", "), "\n")
    }
    cat("Using value column:", value_col[1], "\n")

    # Build expression matrix
    # Use first GSM to get gene IDs
    id_col <- colnames(dt)[1]  # Usually ID_REF or similar
    gene_ids <- dt[[id_col]]

    expr_248423 <- matrix(NA, nrow = length(gene_ids), ncol = length(gsm_list))
    rownames(expr_248423) <- gene_ids
    colnames(expr_248423) <- all_gsm_names

    cat("Building expression matrix for", length(gsm_list), "samples...\n")
    pb <- txtProgressBar(min = 0, max = length(gsm_list), style = 3)
    for (i in seq_along(gsm_list)) {
      gsm_dt <- Table(gsm_list[[i]])
      if (nrow(gsm_dt) > 0 && value_col[1] %in% colnames(gsm_dt)) {
        expr_248423[, i] <- as.numeric(gsm_dt[[value_col[1]]])
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)

    cat("\nExpression matrix:", nrow(expr_248423), "x", ncol(expr_248423), "\n")
    cat("Non-NA values:", sum(!is.na(expr_248423)), "\n")

    # Get metadata
    meta_248423 <- read.csv("results/meta_GSE248423.csv", stringsAsFactors = FALSE)
    cat("Metadata:", nrow(meta_248423), "samples\n")
    if ("disease.state.ch1" %in% colnames(meta_248423)) {
      cat("Disease state distribution:\n")
      print(table(meta_248423$disease.state.ch1))
    }

    save(expr_248423, meta_248423, file = "results/GSE248423_processed.RData")
    cat("Saved: results/GSE248423_processed.RData\n")

  } else {
    cat("No data table in GSM entries — expression stored elsewhere.\n")
    cat("This dataset's counts may only be available as individual sample files.\n")
    cat("Skipping GSE248423 for now.\n")
  }
}


###############################################################################
# SECTION 4: Overall Status + What to Download in Safari
###############################################################################

cat("\n\n")
cat("================================================================\n")
cat("  CURRENT STATUS\n")
cat("================================================================\n\n")

cat("READY TO ANALYZE:\n")
cat("  GSE140829: 587 samples (Illumina BeadChip, Japan)\n")
cat("    Diagnosis: AD=204, Control=249, MCI=134\n")
cat("    Rich metadata: age, sex, APOE, batch, connectivity score\n\n")

cat("  GSE270454: 45 samples (RNA-seq)\n")
cat("    Conditions: AD=10, ASM=11, ASO=14, MCI=10\n")
cat("    Note: No healthy controls — ASM/ASO may serve as comparison\n\n")

if (exists("expr_248423") && sum(!is.na(expr_248423)) > 0) {
  cat("  GSE248423: 196 samples (RNA-seq, extracted from SOFT)\n")
  cat("    Conditions: AD and Control\n\n")
}

cat("DOWNLOAD IN SAFARI (paste these URLs in your browser):\n\n")

cat("  1. GSE63060 (~60 MB):\n")
cat("     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/GSE63060_series_matrix.txt.gz\n")
cat("     Save to: data/raw/GSE63060_series_matrix.txt.gz\n\n")

cat("  2. GSE63061 (~60 MB):\n")
cat("     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/GSE63061_series_matrix.txt.gz\n")
cat("     Save to: data/raw/GSE63061_series_matrix.txt.gz\n\n")

cat("  3. GSE140829 FULL series matrix (~216 MB) — for complete metadata:\n")
cat("     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140829/matrix/GSE140829_series_matrix.txt.gz\n")
cat("     Save to: data/raw/GSE140829_series_matrix.txt.gz\n\n")

cat("After downloading, run this in R to load them:\n")
cat('  gse63060 <- getGEO(filename = "data/raw/GSE63060_series_matrix.txt.gz")\n')
cat('  gse63061 <- getGEO(filename = "data/raw/GSE63061_series_matrix.txt.gz")\n\n')

cat("================================================================\n")
cat("  MINIMUM VIABLE PAPER: GSE140829 alone (587 samples) is enough\n")
cat("  to demonstrate RS/FM framework. Additional datasets strengthen it.\n")
cat("================================================================\n")
