#!/bin/bash

## VIASH START
## VIASH END

# source the helpers & setup env
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

# --- create and validate test-data ---
log "Start creating & fetching test data"

test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"/{fastq_reads,output_test_{1,2,3,4,5}}

# create and validate the fastq reads
create_test_fastq "$test_data_dir/fastq_reads/read_1.fastq" 1 100
check_file_exists "$test_data_dir/fastq_reads/read_1.fastq" "test FASTQ read_1 file"

create_test_fastq "$test_data_dir/fastq_reads/read_2.fastq" 1 100
check_file_exists "$test_data_dir/fastq_reads/read_2.fastq" "test FASTQ read_2 file"

# opted here to fetch data from the HISAT2 repo because of the specific file format
# prefer it to not be dependent on the hisat2-build command for isolation
wget -qO- https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.2.tar.gz \
    | tar -xz -C "$test_data_dir" \
    --strip-components=2 \
    hisat2-2.2.2/example/index \
    hisat2-2.2.2/example/reads

# check if fetched data exists
for i in $(seq 1 8); do
  check_file_exists "$test_data_dir/index/22_20-21M_snp.${i}.ht2" "HISAT2 index file ${i}"
done
check_file_exists "$test_data_dir/reads/reads_1.fa" "reads file 1"
check_file_exists "$test_data_dir/reads/reads_2.fa" "reads file 2"

log "Test data created, fetched and validated"

# --- TEST_1: basic functionality ---
log "Starting TEST_1: basic functionality, running in single-end mode"
log "Now executing $meta_name with basic parameters..."
"$meta_executable" \
  --index_dir "$test_data_dir/index" \
  --index_prefix 22_20-21M_snp \
  --input "$test_data_dir/fastq_reads/read_1.fastq" \
  --output_sam "$test_data_dir/output_test_1/output.sam"

check_file_exists "$test_data_dir/output_test_1/output.sam" "Output SAM file"
check_file_not_empty "$test_data_dir/output_test_1/output.sam" "Output SAM file"
check_file_matches_regex "$test_data_dir/output_test_1/output.sam" "^@(HD|SQ|PG)" "TEST_1 SAM headers"

log "✓ TEST_1 completed successfully"

# --- TEST_2: testing basic functionality in paired-end mode ---
log "Starting TEST_2: basic functionality, running in paired-end mode"
log "Now executing $meta_name with basic parameters..."
"$meta_executable" \
  --index_dir "$test_data_dir/index" \
  --index_prefix 22_20-21M_snp \
  --input "$test_data_dir/fastq_reads/read_1.fastq" \
  --input_r2 "$test_data_dir/fastq_reads/read_2.fastq" \
  --output_sam "$test_data_dir/output_test_2/output.sam"

check_file_exists "$test_data_dir/output_test_2/output.sam" "Output SAM file"
check_file_not_empty "$test_data_dir/output_test_2/output.sam" "Output SAM file"
check_file_matches_regex "$test_data_dir/output_test_2/output.sam" "^@(HD|SQ|PG)" "TEST_2 SAM headers"

log "✓ TEST_2 completed successfully"

# --- TEST_3: testing fasta flag ---
log "Starting TEST_3: testing fasta flag"
log "Now executing $meta_name with basic parameters..."
"$meta_executable" \
  --index_dir "$test_data_dir/index" \
  --index_prefix 22_20-21M_snp \
  --fasta \
  --input "$test_data_dir/reads/reads_1.fa" \
  --input_r2 "$test_data_dir/reads/reads_2.fa" \
  --output_sam "$test_data_dir/output_test_3/output.sam"

check_file_exists "$test_data_dir/output_test_3/output.sam" "Output SAM file"
check_file_not_empty "$test_data_dir/output_test_3/output.sam" "Output SAM file"
check_file_matches_regex "$test_data_dir/output_test_3/output.sam" "^@(HD|SQ|PG)" "TEST_3 SAM headers"

log "✓ TEST_3 completed successfully"

# --- TEST_4: testing --summary_file output ---
log "Starting TEST_4: testing --summary_file optional output"
"$meta_executable" \
  --index_dir "$test_data_dir/index" \
  --index_prefix 22_20-21M_snp \
  --input "$test_data_dir/fastq_reads/read_1.fastq" \
  --output_sam "$test_data_dir/output_test_4/output.sam" \
  --summary_file "$test_data_dir/output_test_4/summary.txt"

check_file_exists "$test_data_dir/output_test_4/output.sam" "Output SAM file"
check_file_exists "$test_data_dir/output_test_4/summary.txt" "Alignment summary file"
check_file_not_empty "$test_data_dir/output_test_4/summary.txt" "Alignment summary file"
check_file_matches_regex "$test_data_dir/output_test_4/summary.txt" "reads; of these" "TEST_4 summary alignment stats"

log "✓ TEST_4 completed successfully"

# --- TEST_5: testing --mp/--sp/--rdg/--rfg multiple-value scoring options ---
log "Starting TEST_5: testing --mp/--sp/--rdg/--rfg (multiple: true integer pairs)"
"$meta_executable" \
  --index_dir "$test_data_dir/index" \
  --index_prefix 22_20-21M_snp \
  --input "$test_data_dir/fastq_reads/read_1.fastq" \
  --output_sam "$test_data_dir/output_test_5/output.sam" \
  --mp 6 --mp 2 \
  --sp 1 --sp 0 \
  --rdg 5 --rdg 3 \
  --rfg 5 --rfg 3

check_file_exists "$test_data_dir/output_test_5/output.sam" "Output SAM file"
check_file_not_empty "$test_data_dir/output_test_5/output.sam" "Output SAM file"
# Confirm Viash's ";"-joined values were correctly converted to hisat2's expected "," format
# (avoid leading "--" in the search pattern so grep doesn't parse it as an option)
check_file_contains "$test_data_dir/output_test_5/output.sam" "mp 6,2" "TEST_5 --mp passed as 6,2"
check_file_contains "$test_data_dir/output_test_5/output.sam" "sp 1,0" "TEST_5 --sp passed as 1,0"
check_file_contains "$test_data_dir/output_test_5/output.sam" "rdg 5,3" "TEST_5 --rdg passed as 5,3"
check_file_contains "$test_data_dir/output_test_5/output.sam" "rfg 5,3" "TEST_5 --rfg passed as 5,3"

log "✓ TEST_5 completed successfully"

# --- TEST_6: a non-integer --mp value must be rejected before hisat2 even runs ---
log "Starting TEST_6: testing --mp rejects a non-integer value"
set +e
"$meta_executable" \
  --index_dir "$test_data_dir/index" \
  --index_prefix 22_20-21M_snp \
  --input "$test_data_dir/fastq_reads/read_1.fastq" \
  --output_sam "$test_data_dir/output_test_6/output.sam" \
  --mp not_an_integer --mp 2 2>"$test_data_dir/output_test_6.err"
exit_code=$?
set -e

if [ "$exit_code" -eq 0 ]; then
  log_error "TEST_6 FAILED: $meta_name should have rejected a non-integer --mp value but succeeded"
  exit 1
fi
check_file_matches_regex "$test_data_dir/output_test_6.err" "mp.*has to be an integer" "TEST_6 error message"

log "✓ TEST_6 completed successfully"

print_test_summary "$meta_name tests"