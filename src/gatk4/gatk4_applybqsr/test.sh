#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the gatk4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

#############################################
# Shared fixtures: reference genome + BAM + BQSR recal table
#############################################

create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

sam_file="$test_data_dir/reads.sam"
create_test_sam_reads "$sam_file" seq1 2000 rg1 sample1 20 100
check_file_exists "$sam_file" "test SAM reads"

sort_and_index_bam "$sam_file" "$test_data_dir/sample.sorted.bam"
check_file_exists "$test_data_dir/sample.sorted.bam" "sorted test BAM"
check_file_exists "$test_data_dir/sample.sorted.bai" "sorted test BAM index"

known_sites_vcf="$test_data_dir/known_sites.vcf"
create_test_known_sites_vcf "$known_sites_vcf" seq1 2000 500
check_file_exists "$known_sites_vcf" "known-sites VCF"

log "Running BaseRecalibrator to produce a recal table..."
recal_table="$test_data_dir/recal.table"
gatk BaseRecalibrator \
  --input "$test_data_dir/sample.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --known-sites "$known_sites_vcf" \
  --output "$recal_table" \
  --verbosity ERROR
check_file_exists "$recal_table" "BQSR recalibration table"
check_file_not_empty "$recal_table" "BQSR recalibration table"

# --- Test Case 1: Basic ApplyBQSR ---
log "Starting TEST 1: Basic ApplyBQSR"

mkdir -p "$meta_temp_dir/test1" "$meta_temp_dir/test1_index"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test1/output.bam" \
  --output_index "$meta_temp_dir/test1_index/recalibrated.bai"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/test1/output.bam" "recalibrated BAM"
check_file_not_empty "$meta_temp_dir/test1/output.bam" "recalibrated BAM"
check_file_exists "$meta_temp_dir/test1_index/recalibrated.bai" "output BAM index"
check_file_not_empty "$meta_temp_dir/test1_index/recalibrated.bai" "output BAM index"
check_file_not_exists "$meta_temp_dir/test1/output.bai" "index next to the output BAM (moved to --output_index)"

log "Validating recalibrated BAM with ValidateSamFile..."
gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test1/output.bam" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test1/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test1/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: ApplyBQSR with options ---
log "Starting TEST 2: ApplyBQSR with intervals, static quantized quals and preserve-qscores-less-than"

mkdir -p "$meta_temp_dir/test2"

intervals_bed="$test_data_dir/intervals.bed"
echo -e "seq1\t0\t1000" > "$intervals_bed"
check_file_exists "$intervals_bed" "intervals BED file"

log "Executing $meta_name with advanced parameters..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test2/output.bam" \
  --output_index "$meta_temp_dir/test2/output.bai" \
  --intervals "$intervals_bed" \
  --static_quantized_quals 10 \
  --static_quantized_quals 20 \
  --static_quantized_quals 30 \
  --preserve_qscores_less_than 20

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/test2/output.bam" "recalibrated BAM with options"
check_file_not_empty "$meta_temp_dir/test2/output.bam" "recalibrated BAM with options"
check_file_exists "$meta_temp_dir/test2/output.bai" "output BAM index (same path GATK writes)"

log "Validating recalibrated BAM with ValidateSamFile..."
gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test2/output.bam" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test2/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test2/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: ApplyBQSR without a reference ---
log "Starting TEST 3: ApplyBQSR without --reference"

mkdir -p "$meta_temp_dir/test3"

log "Executing $meta_name without a reference..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test3/output.bam" \
  --output_index "$meta_temp_dir/test3/output.bai"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/test3/output.bam" "recalibrated BAM without reference"
check_file_not_empty "$meta_temp_dir/test3/output.bam" "recalibrated BAM without reference"
check_file_exists "$meta_temp_dir/test3/output.bai" "output BAM index without reference"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Shared GATK engine options ---
log "Starting TEST 4: Shared GATK engine options"

mkdir -p "$meta_temp_dir/test4"

log "Executing $meta_name with --interval_padding, --sites_only_vcf_output, --create_output_variant_index false, --create_output_bam_index false..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test4/output.bam" \
  --interval_padding 10 \
  --sites_only_vcf_output \
  --create_output_variant_index false \
  --create_output_bam_index false

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/test4/output.bam" "recalibrated BAM (engine options)"
check_file_not_empty "$meta_temp_dir/test4/output.bam" "recalibrated BAM (engine options)"
check_file_not_exists "$meta_temp_dir/test4/output.bai" "output BAM index (should not exist with --create_output_bam_index false)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Missing --output_index ---
log "Starting TEST 5: --output_index is required unless --create_output_bam_index is false"

mkdir -p "$meta_temp_dir/test5"

if "$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test5/output.bam" > "$meta_temp_dir/test5/run.log" 2>&1; then
  log_error "✗ $meta_name did not fail without --output_index"
  exit 1
fi
check_file_contains "$meta_temp_dir/test5/run.log" "Error: --output_index is required unless --create_output_bam_index is false." "error message"
check_file_not_exists "$meta_temp_dir/test5/output.bam" "output BAM (should not be written)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Multiple interval files with exclusions ---
log "Starting TEST 6: Multiple --intervals files (BED and bgzipped VCF) with --exclude_intervals"

mkdir -p "$meta_temp_dir/test6"

# Reads start at 1, 101, ..., 1901 (100bp each). The area is the union
# of seq1:1-600 (BED) and seq1:1550 (VCF), minus seq1:201-300 (BED).
printf 'seq1\t0\t600\n' > "$test_data_dir/intervals_a.bed"
printf 'seq1\t200\t300\n' > "$test_data_dir/exclude.bed"
create_test_known_sites_vcf "$test_data_dir/intervals_site.vcf" seq1 2000 1550
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

log "Executing $meta_name with two --intervals files, --intervals_index and --exclude_intervals..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test6/output.bam" \
  --output_index "$meta_temp_dir/test6/output.bai" \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule UNION \
  --exclude_intervals "$test_data_dir/exclude.bed"

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/test6/output.bam" "recalibrated BAM (multiple intervals)"
check_file_exists "$meta_temp_dir/test6/output.bai" "output BAM index (multiple intervals)"

# ApplyBQSR writes the reads that overlap the area
read_starts=$(samtools view "$meta_temp_dir/test6/output.bam" | cut -f4 | tr '\n' ' ')
expected_starts="1 101 301 401 501 1501 "
if [[ "$read_starts" == "$expected_starts" ]]; then
  log "✓ output BAM only has the reads in the selected territory: $read_starts"
else
  log_error "✗ output BAM read starts are '$read_starts', expected '$expected_starts'"
  exit 1
fi

log "✅ TEST 6 completed successfully"

# --- Test Case 7: CRAM input with a separate index ---
log "Starting TEST 7: CRAM --input with its .crai in another directory"

mkdir -p "$meta_temp_dir/test7"

create_test_cram "$test_data_dir/sample.sorted.bam" "$test_data_dir/reference.fasta" \
  "$test_data_dir/cram/sample.cram" "$test_data_dir/cram_index/sample.cram.crai"
check_file_not_exists "$test_data_dir/cram/sample.cram.crai" "index next to the CRAM input"

log "Executing $meta_name with a CRAM --input and --read_index..."
"$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample.cram.crai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test7/output.bam" \
  --output_index "$meta_temp_dir/test7/output.bai"

log "Validating TEST 7 outputs..."
check_file_exists "$meta_temp_dir/test7/output.bam" "recalibrated BAM from CRAM input"
check_file_exists "$meta_temp_dir/test7/output.bai" "output BAM index from CRAM input"

gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test7/output.bam" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test7/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test7/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

read_count=$(samtools view -c "$meta_temp_dir/test7/output.bam")
if [[ "$read_count" == "20" ]]; then
  log "✓ output BAM has all $read_count reads of the CRAM input"
else
  log_error "✗ output BAM has $read_count reads, expected 20"
  exit 1
fi

log "✅ TEST 7 completed successfully"

# --- Test Case 8: CRAM output ---
log "Starting TEST 8: CRAM --output with --output_index in another directory"

mkdir -p "$meta_temp_dir/test8" "$meta_temp_dir/test8_index"

log "Executing $meta_name with a CRAM --output..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test8/output.cram" \
  --output_index "$meta_temp_dir/test8_index/output.cram.bai"

log "Validating TEST 8 outputs..."
check_file_exists "$meta_temp_dir/test8/output.cram" "recalibrated CRAM"
check_file_exists "$meta_temp_dir/test8_index/output.cram.bai" "output CRAM index"
check_file_not_exists "$meta_temp_dir/test8/output.cram.bai" "index next to the output CRAM (moved to --output_index)"
check_file_not_exists "$meta_temp_dir/test8/output.cram.crai" "CRAI index next to the output CRAM"
check_file_not_exists "$meta_temp_dir/test8/output.bai" "BAI index next to the output CRAM"

if [[ "$(head -c 4 "$meta_temp_dir/test8/output.cram")" == "CRAM" ]]; then
  log "✓ output has the CRAM magic bytes"
else
  log_error "✗ output does not have the CRAM magic bytes"
  exit 1
fi
if [[ "$(head -c 4 "$meta_temp_dir/test8_index/output.cram.bai")" == $'BAI\x01' ]]; then
  log "✓ output index has the BAI magic bytes"
else
  log_error "✗ output index does not have the BAI magic bytes"
  exit 1
fi

gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test8/output.cram" \
  --REFERENCE_SEQUENCE "$test_data_dir/reference.fasta" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test8/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test8/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: CRAM input without a reference ---
log "Starting TEST 9: CRAM --input without --reference"

mkdir -p "$meta_temp_dir/test9"

if "$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample.cram.crai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test9/output.bam" \
  --output_index "$meta_temp_dir/test9/output.bai" > "$meta_temp_dir/test9/run.log" 2>&1; then
  log_error "✗ $meta_name did not fail with a CRAM --input and no --reference"
  exit 1
fi
check_file_contains "$meta_temp_dir/test9/run.log" "Error: --reference is required for CRAM input or output." "error message"
check_file_not_exists "$meta_temp_dir/test9/output.bam" "output BAM (should not be written)"

log "✅ TEST 9 completed successfully"

print_test_summary "All tests"
