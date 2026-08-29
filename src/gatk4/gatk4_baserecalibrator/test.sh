#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the GATK4-specific test helpers
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
# Shared fixtures: reference genome + BAM
#############################################

create_test_reference "$test_data_dir/ref.fasta" 1 1000
check_file_exists "$test_data_dir/ref.fasta" "test reference genome"
check_file_exists "$test_data_dir/ref.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/ref.fasta.fai" "reference FASTA index"

sam_file="$test_data_dir/reads.sam"
create_test_sam_reads "$sam_file" seq1 1000 rg1 sample1 10 100
check_file_exists "$sam_file" "test SAM file"

sort_and_index_bam "$sam_file" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bam" "sorted BAM file"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"

# --- Test Case 1: Basic functionality with a single known-sites file ---
log "Starting TEST 1: Basic functionality with a single known-sites file"

known_sites1="$test_data_dir/known_sites1.vcf"
create_test_known_sites_vcf "$known_sites1" seq1 1000 50
check_file_exists "$known_sites1" "known-sites VCF file"
check_file_exists "$known_sites1.idx" "known-sites VCF index file"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal1.table"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/recal1.table" "recalibration table"
check_file_not_empty "$meta_temp_dir/recal1.table" "recalibration table"
check_file_contains "$meta_temp_dir/recal1.table" "#:GATKReport" "recalibration table GATKReport header"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Multiple known-sites files and custom options ---
log "Starting TEST 2: Multiple known-sites files and custom options"

known_sites2="$test_data_dir/known_sites2.vcf"
create_test_known_sites_vcf "$known_sites2" seq1 1000 750
check_file_exists "$known_sites2" "second known-sites VCF file"
check_file_exists "$known_sites2.idx" "second known-sites VCF index file"

log "Executing $meta_name with multiple known-sites files and custom options..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$known_sites1" \
  --known_sites "$known_sites2" \
  --low_quality_tail 5 \
  --quantizing_levels 10 \
  --indels_context_size 2 \
  --mismatches_context_size 1 \
  --default_base_qualities 30 \
  --mismatches_default_quality 40 \
  --preserve_qscores_less_than 10 \
  --use_original_qualities \
  --output "$meta_temp_dir/recal2.table"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/recal2.table" "recalibration table with custom options"
check_file_not_empty "$meta_temp_dir/recal2.table" "recalibration table with custom options"
check_file_contains "$meta_temp_dir/recal2.table" "#:GATKReport" "recalibration table GATKReport header"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Shared GATK engine options ---
log "Starting TEST 3: Shared GATK engine options"

log "Executing $meta_name with --interval_padding, --read_filter and --disable_tool_default_read_filters..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$known_sites1" \
  --interval_padding 10 \
  --read_filter MappingQualityReadFilter \
  --disable_tool_default_read_filters \
  --output "$meta_temp_dir/recal_engine.table"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/recal_engine.table" "recalibration table (engine options)"
check_file_not_empty "$meta_temp_dir/recal_engine.table" "recalibration table (engine options)"
check_file_contains "$meta_temp_dir/recal_engine.table" "#:GATKReport" "recalibration table GATKReport header (engine options)"

log "✅ TEST 3 completed successfully"

print_test_summary "All tests"
