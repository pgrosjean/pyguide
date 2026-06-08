#!/usr/bin/env bash
# Download GuideScan2 hg38 index (~2.2 GB) and hg38 chrom sizes.
# Usage: bash scripts/download_guidescan_data.sh [--data-dir /path/to/dir]

set -euo pipefail

DATA_DIR="${HOME}/.pyguide_data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Usage: $0 [--data-dir /path/to/dir]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${DATA_DIR}"

INDEX_URL="https://guidescan.com/indices/hg38.zip"
INDEX_ZIP="${DATA_DIR}/hg38.zip"

CHROM_SIZES_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
CHROM_SIZES_DEST="$(dirname "$0")/../data/hg38.chrom.sizes"

echo "=== GuideScan2 hg38 Data Download ==="
echo "Data directory: ${DATA_DIR}"
echo ""

# Download hg38 GuideScan2 index
echo "[DOWNLOAD] hg38 GuideScan2 index (~2.2 GB) ..."
echo "  Note: Using --insecure due to expired SSL cert on guidescan.com"
curl --insecure -L --progress-bar -o "${INDEX_ZIP}" "${INDEX_URL}"
echo "[EXTRACT] Extracting hg38.zip ..."
unzip -o "${INDEX_ZIP}" -d "${DATA_DIR}"
EXTRACTED=$(find "${DATA_DIR}" -not -name "*.zip" -newer "${INDEX_ZIP}" -type f 2>/dev/null | head -1)
if [[ -z "${EXTRACTED}" ]]; then
    echo "ERROR: Could not find index file after extraction." >&2
    echo "Contents of ${DATA_DIR}:" >&2
    ls "${DATA_DIR}" >&2
    exit 1
fi
echo "[OK] Index extracted: ${EXTRACTED}"

# Download hg38 chrom sizes (for BigWig output)
echo ""
DEST_DIR="$(dirname "${CHROM_SIZES_DEST}")"
mkdir -p "${DEST_DIR}"
if [[ -f "${CHROM_SIZES_DEST}" ]]; then
    echo "[SKIP] hg38.chrom.sizes already exists: ${CHROM_SIZES_DEST}"
else
    echo "[DOWNLOAD] hg38 chromosome sizes ..."
    curl -L --progress-bar -o "${CHROM_SIZES_DEST}" "${CHROM_SIZES_URL}"
    echo "[OK] chrom sizes written: ${CHROM_SIZES_DEST}"
fi

echo ""
echo "=== Done ==="
echo "Index location: ${DATA_DIR}"
echo "Run 'ls ${DATA_DIR}' to see the extracted index filename."
echo "Pass to pyguide-tiling with:  --index ${DATA_DIR}/<index-filename>"
