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
  --read_index "$test_data_dir/reads.sorted.bai" \
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
  --read_index "$test_data_dir/reads.sorted.bai" \
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
  --read_index "$test_data_dir/reads.sorted.bai" \
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

# --- Test Case 4: Bgzipped known-sites files with separate indexes ---
log "Starting TEST 4: Bgzipped known-sites files with separate indexes"

create_test_vcf_gz_with_separate_index "$known_sites1" \
  "$test_data_dir/known_sites1.vcf.gz" "$test_data_dir/indexes/known_sites1.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$known_sites2" \
  "$test_data_dir/known_sites2.vcf.gz" "$test_data_dir/indexes/known_sites2.vcf.gz.tbi"

log "Executing $meta_name with --known_sites_index..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$test_data_dir/known_sites1.vcf.gz" \
  --known_sites "$test_data_dir/known_sites2.vcf.gz" \
  --known_sites_index "$test_data_dir/indexes/known_sites1.vcf.gz.tbi" \
  --known_sites_index "$test_data_dir/indexes/known_sites2.vcf.gz.tbi" \
  --output "$meta_temp_dir/recal_gz.table"

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/recal_gz.table" "recalibration table (bgzipped known sites)"
check_file_not_empty "$meta_temp_dir/recal_gz.table" "recalibration table (bgzipped known sites)"
check_file_contains "$meta_temp_dir/recal_gz.table" "#:GATKReport" "recalibration table GATKReport header (bgzipped known sites)"
check_file_not_exists "$test_data_dir/known_sites1.vcf.gz.tbi" "index next to the first known-sites input"
check_file_not_exists "$test_data_dir/known_sites2.vcf.gz.tbi" "index next to the second known-sites input"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Bgzipped known-sites file without an index ---
log "Starting TEST 5: Bgzipped known-sites file without an index"

log "Executing $meta_name with a .vcf.gz known-sites file and no --known_sites_index..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$test_data_dir/known_sites1.vcf.gz" \
  --output "$meta_temp_dir/recal_no_index.table" 2> "$meta_temp_dir/test5.log"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/recal_no_index.table" "recalibration table (no known-sites index)"
check_file_contains "$meta_temp_dir/recal_no_index.table" "#:GATKReport" "recalibration table GATKReport header (no known-sites index)"
check_file_contains "$meta_temp_dir/test5.log" "no index was provided" "warning for the missing index"
check_file_not_exists "$test_data_dir/known_sites1.vcf.gz.tbi" "index next to the known-sites input"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Fewer --known_sites_index files than --known_sites files ---
log "Starting TEST 6: Fewer --known_sites_index files than --known_sites files"

log "Executing $meta_name with two known-sites files and one index (should fail)..."
if "$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$test_data_dir/known_sites1.vcf.gz" \
  --known_sites "$test_data_dir/known_sites2.vcf.gz" \
  --known_sites_index "$test_data_dir/indexes/known_sites1.vcf.gz.tbi" \
  --output "$meta_temp_dir/recal_bad.table" > "$meta_temp_dir/test6.log" 2>&1; then
  log_error "✗ $meta_name did not fail with fewer --known_sites_index files than --known_sites files"
  exit 1
fi
check_file_contains "$meta_temp_dir/test6.log" "Error: --known_sites_index must be given once for each file in --known_sites" "error message"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Intersection of multiple interval files ---
log "Starting TEST 7: Multiple --intervals files (BED and bgzipped VCF) with --interval_set_rule INTERSECTION"

# Reads start at 1, 101, ..., 901 (100bp each). The area is the
# intersection of seq1:401-800 (BED) and seq1:650 (VCF), so only the read at
# 601 is used (100 bases, no known site in it). With UNION, the reads at 401,
# 501, 601 and 701 would be used.
printf 'seq1\t400\t800\n' > "$test_data_dir/intervals_a.bed"
create_test_known_sites_vcf "$test_data_dir/intervals_site.vcf" seq1 1000 650
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

log "Executing $meta_name with two --intervals files, --intervals_index and --interval_set_rule INTERSECTION..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal_intervals.table" \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule INTERSECTION

log "Validating TEST 7 outputs..."
check_file_exists "$meta_temp_dir/recal_intervals.table" "recalibration table (multiple intervals)"

# The read group row of RecalTable0 has the number of observed bases
observations=$(awk '/^#:GATKTable:RecalTable0:/ {in_table = 1; next} in_table && $1 == "rg1" {print $5; exit}' \
  "$meta_temp_dir/recal_intervals.table")
if [[ "$observations" == "100" ]]; then
  log "✓ recalibration table only has the bases of the read in the selected territory ($observations observations)"
else
  log_error "✗ recalibration table has '$observations' observations, expected 100"
  exit 1
fi

log "✅ TEST 7 completed successfully"

# --- Test Case 8: CRAM input with a separate index ---
log "Starting TEST 8: CRAM --input with its .crai in another directory"

create_test_cram "$test_data_dir/reads.sorted.bam" "$test_data_dir/ref.fasta" \
  "$test_data_dir/cram/reads.cram" "$test_data_dir/cram_index/reads.cram.crai"
check_file_not_exists "$test_data_dir/cram/reads.cram.crai" "index next to the CRAM input"

log "Executing $meta_name with a CRAM --input and --read_index..."
"$meta_executable" \
  --input "$test_data_dir/cram/reads.cram" \
  --read_index "$test_data_dir/cram_index/reads.cram.crai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal_cram.table"

log "Validating TEST 8 outputs..."
check_file_exists "$meta_temp_dir/recal_cram.table" "recalibration table (CRAM input)"
check_file_contains "$meta_temp_dir/recal_cram.table" "#:GATKReport" "recalibration table GATKReport header (CRAM input)"

# The CRAM has the same reads as the BAM, so the read group row of
# RecalTable0 must be the same as for the BAM input in TEST 1
observations() {
  awk '/^#:GATKTable:RecalTable0:/ {in_table = 1; next} in_table && $1 == "rg1" {print $5; exit}' "$1"
}
cram_observations=$(observations "$meta_temp_dir/recal_cram.table")
bam_observations=$(observations "$meta_temp_dir/recal1.table")
if [[ -n "$cram_observations" && "$cram_observations" == "$bam_observations" ]]; then
  log "✓ recalibration table from the CRAM input has the same observations as from the BAM input ($cram_observations)"
else
  log_error "✗ recalibration table from the CRAM input has '$cram_observations' observations, expected '$bam_observations'"
  exit 1
fi

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Two input files ---
log "Starting TEST 9: Two --input files give one recalibration table"

# A second BAM with another read group of the same sample, as for a second
# lane (5 reads at 51, 251, 451, 651 and 851)
create_test_sam "$test_data_dir/reads_rg2.sam" seq1 1000 rg2 sample1 51 251 451 651 851
sed -i 's/^read/rg2_read/' "$test_data_dir/reads_rg2.sam"
sort_and_index_bam "$test_data_dir/reads_rg2.sam" "$test_data_dir/reads_rg2.sorted.bam"
check_file_exists "$test_data_dir/reads_rg2.sorted.bai" "second BAM index"

log "Executing $meta_name with two --input and --read_index files..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  --input "$test_data_dir/reads_rg2.sorted.bam" \
  --read_index "$test_data_dir/reads_rg2.sorted.bai" \
  --reference "$test_data_dir/ref.fasta" \
  --reference_fai "$test_data_dir/ref.fasta.fai" \
  --reference_dict "$test_data_dir/ref.dict" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal_two_inputs.table"

log "Validating TEST 9 outputs..."
check_file_contains "$meta_temp_dir/recal_two_inputs.table" "#:GATKReport" "recalibration table GATKReport header (two inputs)"

# RecalTable0 has one row for each read group. The rg1 row must be the same
# as for the first BAM alone (TEST 1), and the rg2 row has the bases of the 5
# reads of the second BAM (5 x 100, none of them has the known site).
rg_observations() {
  awk -v rg="$2" '/^#:GATKTable:RecalTable0:/ {in_table = 1; next} in_table && $1 == rg {print $5; exit}' "$1"
}
rg1_observations=$(rg_observations "$meta_temp_dir/recal_two_inputs.table" rg1)
rg2_observations=$(rg_observations "$meta_temp_dir/recal_two_inputs.table" rg2)
expected_rg1_observations=$(rg_observations "$meta_temp_dir/recal1.table" rg1)
if [[ -n "$rg1_observations" && "$rg1_observations" == "$expected_rg1_observations" ]]; then
  log "✓ recalibration table has the rg1 observations of the first BAM ($rg1_observations)"
else
  log_error "✗ recalibration table has '$rg1_observations' rg1 observations, expected '$expected_rg1_observations'"
  exit 1
fi
if [[ "$rg2_observations" == "500" ]]; then
  log "✓ recalibration table has the rg2 observations of the second BAM ($rg2_observations)"
else
  log_error "✗ recalibration table has '$rg2_observations' rg2 observations, expected 500"
  exit 1
fi

log "✅ TEST 9 completed successfully"

print_test_summary "All tests"
