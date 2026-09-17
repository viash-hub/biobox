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

mkdir -p "$meta_temp_dir/test1"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --bai "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test1/output.bam"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/test1/output.bam" "recalibrated BAM"
check_file_not_empty "$meta_temp_dir/test1/output.bam" "recalibrated BAM"

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
  --bai "$test_data_dir/sample.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test2/output.bam" \
  --intervals "$intervals_bed" \
  --static_quantized_quals 10 \
  --static_quantized_quals 20 \
  --static_quantized_quals 30 \
  --preserve_qscores_less_than 20

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/test2/output.bam" "recalibrated BAM with options"
check_file_not_empty "$meta_temp_dir/test2/output.bam" "recalibrated BAM with options"

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
  --bai "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test3/output.bam"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/test3/output.bam" "recalibrated BAM without reference"
check_file_not_empty "$meta_temp_dir/test3/output.bam" "recalibrated BAM without reference"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Shared GATK engine options ---
log "Starting TEST 4: Shared GATK engine options"

mkdir -p "$meta_temp_dir/test4"

log "Executing $meta_name with --interval_padding, --sites_only_vcf_output, --create_output_variant_index false, --create_output_bam_index false..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --bai "$test_data_dir/sample.sorted.bai" \
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

print_test_summary "All tests"
