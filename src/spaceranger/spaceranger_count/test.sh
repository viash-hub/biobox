#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="target/native/spaceranger/spaceranger_count/spaceranger_count"
meta_resources_dir="src/spaceranger/spaceranger_count"
## VIASH END

test_data="$meta_resources_dir/test_data"

echo "> Default test run"
"$meta_executable" \
    --id test_spaceranger \
    --transcriptome "$test_data/GRCh38" \
    --fastqs "$test_data/subsampled_fastqs" \
    --probe_set "$test_data/Visium_HD_Human_Pancreas_probe_set.csv" \
    --cytaimage "$test_data/Visium_HD_Human_Pancreas_image.tif" \
    --image "$test_data/Visium_HD_Human_Pancreas_tissue_image.btf" \
    --create_bam false \
    --lanes 1 \
    --slide H1-HBNMBMC \
    --area D1 \
    --nosecondary

echo "> Checking outputs..."

# Define output directory
OUT_DIR="test_spaceranger/outs"

# Function to check if file exists and is non-empty
check_file() {
    local file=$1
    local description=$2
    echo -n "Checking $description... "
    if [ ! -f "$file" ]; then
        echo "FAIL (file not found)"
        exit 1
    elif [ ! -s "$file" ]; then
        echo "FAIL (file is empty)"
        exit 1
    else
        echo "OK"
    fi
}

# Function to check if directory exists
check_dir() {
    local dir=$1
    local description=$2
    echo -n "Checking $description... "
    if [ ! -d "$dir" ]; then
        echo "FAIL (directory not found)"
        exit 1
    else
        echo "OK"
    fi
}

# Check summary files
check_file "$OUT_DIR/web_summary.html" "web summary"
check_file "$OUT_DIR/metrics_summary.csv" "metrics summary"

# Check spatial files
check_file "$OUT_DIR/spatial/aligned_fiducials.jpg" "fiducial alignment image"
check_file "$OUT_DIR/spatial/detected_tissue_image.jpg" "tissue detection image"
check_file "$OUT_DIR/spatial/scalefactors_json.json" "scale factors JSON"
check_file "$OUT_DIR/spatial/tissue_positions.csv" "tissue positions"

# Check expression matrices
check_dir "$OUT_DIR/filtered_feature_bc_matrix" "filtered matrix directory"
check_file "$OUT_DIR/filtered_feature_bc_matrix/barcodes.tsv.gz" "filtered barcodes"
check_file "$OUT_DIR/filtered_feature_bc_matrix/features.tsv.gz" "filtered features"
check_file "$OUT_DIR/filtered_feature_bc_matrix/matrix.mtx.gz" "filtered matrix"

check_dir "$OUT_DIR/raw_feature_bc_matrix" "raw matrix directory"
check_file "$OUT_DIR/raw_feature_bc_matrix/barcodes.tsv.gz" "raw barcodes"
check_file "$OUT_DIR/raw_feature_bc_matrix/features.tsv.gz" "raw features"
check_file "$OUT_DIR/raw_feature_bc_matrix/matrix.mtx.gz" "raw matrix"

# Check visualization files
check_file "$OUT_DIR/spatial/spatial.loup" "Loupe file"

# Basic content validation
echo "> Performing basic content validation..."

# Check JSON validity
echo -n "Validating JSON structure... "
if ! jq '.' "$OUT_DIR/spatial/scalefactors_json.json" >/dev/null 2>&1; then
    echo "FAIL (invalid JSON)"
    exit 1
else
    echo "OK"
fi

# Check CSV format of tissue positions
echo -n "Validating tissue positions CSV... "
if ! head -n 1 "$OUT_DIR/spatial/tissue_positions.csv" | grep -q ","; then
    echo "FAIL (invalid CSV format)"
    exit 1
else
    echo "OK"
fi

# Check if matrices are valid gzip files
echo -n "Validating matrix compression... "
if ! gzip -t "$OUT_DIR/filtered_feature_bc_matrix/matrix.mtx.gz" 2>/dev/null; then
    echo "FAIL (invalid gzip for filtered matrix)"
    exit 1
else
    echo "OK"
fi

# Check if web summary contains essential sections
echo -n "Checking web summary content... "
if ! grep -q "Spatial RNA" "$OUT_DIR/web_summary.html"; then
    echo "FAIL (missing essential content)"
    exit 1
else
    echo "OK"
fi

echo "> All tests passed successfully!"