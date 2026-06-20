#!/bin/bash
###############################################################################
# Fix truncated series matrices — Run in Terminal
#
#   cd "."
#   bash fix_downloads_v2.sh
###############################################################################

cd "./data/raw"

echo "=== Deleting truncated files and re-downloading ==="
echo ""

# GSE140829 series matrix: got 2.9 MB, need ~216 MB
echo "[1/3] GSE140829 series matrix (~216 MB — largest, ~5-10 min)..."
rm -f GSE140829_series_matrix.txt.gz
curl -L --retry 5 --retry-delay 15 --connect-timeout 30 --max-time 1800 \
  -o GSE140829_series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140829/matrix/GSE140829_series_matrix.txt.gz"
echo "  Size: $(du -h GSE140829_series_matrix.txt.gz | cut -f1)"
echo ""

# GSE63060 series matrix: got 34.5 MB, need ~60 MB
echo "[2/3] GSE63060 series matrix (~60 MB)..."
rm -f GSE63060_series_matrix.txt.gz
curl -L --retry 5 --retry-delay 15 --connect-timeout 30 --max-time 600 \
  -o GSE63060_series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/GSE63060_series_matrix.txt.gz"
echo "  Size: $(du -h GSE63060_series_matrix.txt.gz | cut -f1)"
echo ""

# GSE63061 series matrix: got 8.7 MB, need ~60 MB
echo "[3/3] GSE63061 series matrix (~60 MB)..."
rm -f GSE63061_series_matrix.txt.gz
curl -L --retry 5 --retry-delay 15 --connect-timeout 30 --max-time 600 \
  -o GSE63061_series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/GSE63061_series_matrix.txt.gz"
echo "  Size: $(du -h GSE63061_series_matrix.txt.gz | cut -f1)"
echo ""

echo "=== Verification ==="
echo "Expected sizes:"
echo "  GSE140829_series_matrix.txt.gz  ~216 MB"
echo "  GSE63060_series_matrix.txt.gz   ~60 MB"
echo "  GSE63061_series_matrix.txt.gz   ~60 MB"
echo ""
echo "Actual sizes:"
ls -lh GSE140829_series_matrix.txt.gz GSE63060_series_matrix.txt.gz GSE63061_series_matrix.txt.gz
echo ""
echo "If still truncated, download in your browser instead:"
echo "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140829/matrix/"
echo "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/"
echo "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/"
echo ""
echo "After download, go back to RStudio and re-run the loading script."
