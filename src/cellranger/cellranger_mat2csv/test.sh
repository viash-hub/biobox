#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting tests for $meta_name"

##############################################################
log "Generating a feature-barcode matrix with cellranger count"
##############################################################
# Cell Ranger ships no example matrix, so create one from the tiny FASTQ files and the tiny
# reference of the Cell Ranger installation. A single run yields both representations that
# mat2csv accepts: the filtered_feature_bc_matrix.h5 file and the filtered_feature_bc_matrix
# MEX folder. The reference holds a single genome, "tiny_ref".
count_dir="$meta_temp_dir/count"
cellranger count \
  --id=run \
  --output-dir="$count_dir" \
  --fastqs=/opt/cellranger-10.0.0/external/cellranger_tiny_fastq \
  --transcriptome=/opt/cellranger-10.0.0/external/cellranger_tiny_ref \
  --create-bam=false \
  --nosecondary \
  --disable-ui \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem="$meta_memory_gb"}

matrix_h5="$count_dir/outs/filtered_feature_bc_matrix.h5"
matrix_mex="$count_dir/outs/filtered_feature_bc_matrix"
check_file_exists "$matrix_h5" "Feature-barcode h5 file"
check_dir_exists "$matrix_mex" "Feature-barcode MEX folder"

# The matrix holds 273 features and 1137 barcodes, so the dense CSV has 274 lines (one
# header plus one line per feature) and 1138 fields per line (one feature id column plus
# one column per barcode).
expected_lines=274
expected_fields=1138

##############################################################
log "Starting TEST 1: MEX folder input"
##############################################################
output_mex="$meta_temp_dir/test1/matrix.csv"
"$meta_executable" \
  --input "$matrix_mex" \
  --output "$output_mex"

check_file_exists "$output_mex" "Dense CSV"
check_file_not_empty "$output_mex" "Dense CSV"
check_file_line_count "$output_mex" "$expected_lines" "Dense CSV"
check_file_matches_regex "$output_mex" "^,[ACGT]+-1," "Dense CSV header"
check_file_contains "$output_mex" "ENSG00000234703" "Dense CSV"
fields=$(head -1 "$output_mex" | awk -F, '{print NF}')
[[ "$fields" -eq "$expected_fields" ]] || {
  log_error "Dense CSV header has $fields fields, expected $expected_fields"
  exit 1
}
log "TEST 1 completed successfully"

##############################################################
log "Starting TEST 2: h5 input"
##############################################################
output_h5="$meta_temp_dir/test2/matrix.csv"
"$meta_executable" \
  --input "$matrix_h5" \
  --output "$output_h5"

check_file_exists "$output_h5" "Dense CSV"
check_file_line_count "$output_h5" "$expected_lines" "Dense CSV"
# Both representations describe the same matrix, so they must convert to the same CSV.
cmp -s "$output_h5" "$output_mex" || {
  log_error "CSV of the h5 input differs from the CSV of the MEX input"
  exit 1
}
log "CSV of the h5 input is identical to the CSV of the MEX input"
log "TEST 2 completed successfully"

##############################################################
log "Starting TEST 3: h5 input, restricted to the genome of the reference"
##############################################################
output_genome="$meta_temp_dir/test3/matrix.csv"
"$meta_executable" \
  --input "$matrix_h5" \
  --output "$output_genome" \
  --genome tiny_ref

check_file_exists "$output_genome" "Dense CSV"
check_file_line_count "$output_genome" "$expected_lines" "Dense CSV"
# Every feature belongs to tiny_ref, the only genome of the reference, so selecting it
# keeps the matrix as it is.
cmp -s "$output_genome" "$output_h5" || {
  log_error "CSV restricted to tiny_ref differs from the unrestricted CSV"
  exit 1
}
log "CSV restricted to tiny_ref is identical to the unrestricted CSV"
log "TEST 3 completed successfully"

##############################################################
log "Starting TEST 4: h5 input, restricted to an absent genome"
##############################################################
output_absent="$meta_temp_dir/test4/matrix.csv"
command_output="$meta_temp_dir/test4_output.txt"
if "$meta_executable" \
  --input "$matrix_h5" \
  --output "$output_absent" \
  --genome GRCh38 > "$command_output" 2>&1; then
  log_error "Expected a non-zero exit code when selecting a genome the matrix does not hold"
  exit 1
fi
log "Selecting a genome the matrix does not hold fails, as expected"
check_file_contains "$command_output" "Genome 'GRCh38' not found" "Command output"
check_file_not_exists "$output_absent" "Dense CSV"
log "TEST 4 completed successfully"

print_test_summary "All tests"
