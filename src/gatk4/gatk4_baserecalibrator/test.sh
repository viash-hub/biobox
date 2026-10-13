#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the GATK4-specific test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

#############################################
# Shared fixtures
#############################################

create_test_reference "$test_data_dir/ref.fasta" 1 1000
check_file_exists "$test_data_dir/ref.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/ref.fasta.fai" "reference FASTA index"
reference_args=(
  --reference "$test_data_dir/ref.fasta"
  --reference_fai "$test_data_dir/ref.fasta.fai"
  --reference_dict "$test_data_dir/ref.dict"
)

# First BAM: 10 reads (rg1) at 1, 101, ..., 901
create_test_sam_reads "$test_data_dir/reads.sam" seq1 1000 rg1 sample1 10 100
sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"

# The same reads as a CRAM, with the .crai in another directory
create_test_cram "$test_data_dir/reads.sorted.bam" "$test_data_dir/ref.fasta" \
  "$test_data_dir/cram/reads.cram" "$test_data_dir/cram_index/reads.cram.crai"
check_file_not_exists "$test_data_dir/cram/reads.cram.crai" "index next to the CRAM input"

# Second BAM: another read group of the same sample, as for a second lane
# (5 reads at 51, 251, 451, 651 and 851)
create_test_sam "$test_data_dir/reads_rg2.sam" seq1 1000 rg2 sample1 51 251 451 651 851
sed -i 's/^read/rg2_read/' "$test_data_dir/reads_rg2.sam"
sort_and_index_bam "$test_data_dir/reads_rg2.sam" "$test_data_dir/reads_rg2.sorted.bam"
check_file_exists "$test_data_dir/reads_rg2.sorted.bai" "second BAM index"

# Known sites: a plain indexed VCF, and bgzipped copies of two VCFs with
# their indexes in another directory
known_sites1="$test_data_dir/known_sites1.vcf"
create_test_known_sites_vcf "$known_sites1" seq1 1000 50
check_file_exists "$known_sites1.idx" "known-sites VCF index file"
create_test_known_sites_vcf --no-index "$test_data_dir/known_sites2.vcf" seq1 1000 750
create_test_vcf_gz_with_separate_index "$known_sites1" \
  "$test_data_dir/known_sites1.vcf.gz" "$test_data_dir/indexes/known_sites1.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/known_sites2.vcf" \
  "$test_data_dir/known_sites2.vcf.gz" "$test_data_dir/indexes/known_sites2.vcf.gz.tbi"

# Print the number of observed bases of one read group (RecalTable0 has one
# row for each read group)
#
# Usage: rg_observations "/path/to/recal.table" read_group
rg_observations() {
  awk -v rg="$2" '/^#:GATKTable:RecalTable0:/ {in_table = 1; next} in_table && $1 == rg {print $5; exit}' "$1"
}

# Check one value of the "Arguments" table of a recalibration report
#
# Usage: check_report_argument "/path/to/recal.table" argument expected_value
check_report_argument() {
  local value
  value=$(awk -v arg="$2" '/^#:GATKTable:Arguments:/ {in_table = 1; next} in_table && /^$/ {exit} in_table && $1 == arg {print $2; exit}' "$1")
  if [[ "$value" == "$3" ]]; then
    log "✓ recalibration report has $2 = $3 in the arguments table"
  else
    log_error "✗ recalibration report has $2 = '$value' in the arguments table, expected '$3': $1"
    exit 1
  fi
}

# --- Test Case 1: Basic run (defaults) ---
log "Starting TEST 1: Basic run (defaults)"

"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  "${reference_args[@]}" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal1.table"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/recal1.table" "recalibration table"
check_file_contains "$meta_temp_dir/recal1.table" "#:GATKReport" "recalibration table GATKReport header"
expected_rg1_observations=$(rg_observations "$meta_temp_dir/recal1.table" rg1)

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Inputs and indexes"

# Run A: a CRAM (separate .crai) and a BAM as two --input files
log "Executing $meta_name with a CRAM and a BAM --input..."
"$meta_executable" \
  --input "$test_data_dir/cram/reads.cram" \
  --read_index "$test_data_dir/cram_index/reads.cram.crai" \
  --input "$test_data_dir/reads_rg2.sorted.bam" \
  --read_index "$test_data_dir/reads_rg2.sorted.bai" \
  "${reference_args[@]}" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal_two_inputs.table"

check_file_contains "$meta_temp_dir/recal_two_inputs.table" "#:GATKReport" "recalibration table GATKReport header (two inputs)"
# The CRAM has the same reads as the BAM of TEST 1, so the rg1 row must be the
# same. The rg2 row has the bases of the 5 reads of the second BAM (5 x 100,
# none of them has the known site).
rg1_observations=$(rg_observations "$meta_temp_dir/recal_two_inputs.table" rg1)
rg2_observations=$(rg_observations "$meta_temp_dir/recal_two_inputs.table" rg2)
if [[ -n "$rg1_observations" && "$rg1_observations" == "$expected_rg1_observations" ]]; then
  log "✓ recalibration table has the rg1 observations of the BAM input of TEST 1 from the CRAM ($rg1_observations)"
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

# Run B: two bgzipped known-sites files with separate indexes, and the
# recalibration options. The options that GATK writes to the report
# arguments table are checked there. --default_base_qualities,
# --preserve_qscores_less_than and --use_original_qualities are not in that
# table (pass-through only).
log "Executing $meta_name with --known_sites_index and recalibration options..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  "${reference_args[@]}" \
  --known_sites "$test_data_dir/known_sites1.vcf.gz" \
  --known_sites "$test_data_dir/known_sites2.vcf.gz" \
  --known_sites_index "$test_data_dir/indexes/known_sites1.vcf.gz.tbi" \
  --known_sites_index "$test_data_dir/indexes/known_sites2.vcf.gz.tbi" \
  --low_quality_tail 5 \
  --quantizing_levels 10 \
  --indels_context_size 2 \
  --mismatches_context_size 1 \
  --default_base_qualities 30 \
  --mismatches_default_quality 40 \
  --preserve_qscores_less_than 10 \
  --use_original_qualities \
  --output "$meta_temp_dir/recal_gz.table"

check_file_not_empty "$meta_temp_dir/recal_gz.table" "recalibration table (bgzipped known sites)"
check_file_contains "$meta_temp_dir/recal_gz.table" "#:GATKReport" "recalibration table GATKReport header (bgzipped known sites)"
check_file_not_exists "$test_data_dir/known_sites1.vcf.gz.tbi" "index next to the first known-sites input"
check_file_not_exists "$test_data_dir/known_sites2.vcf.gz.tbi" "index next to the second known-sites input"
check_report_argument "$meta_temp_dir/recal_gz.table" low_quality_tail 5
check_report_argument "$meta_temp_dir/recal_gz.table" quantizing_levels 10
check_report_argument "$meta_temp_dir/recal_gz.table" indels_context_size 2
check_report_argument "$meta_temp_dir/recal_gz.table" mismatches_context_size 1
check_report_argument "$meta_temp_dir/recal_gz.table" mismatches_default_quality 40

# Run C: a bgzipped known-sites file without an index, and the engine
# options (pass-through)
log "Executing $meta_name with a .vcf.gz known-sites file, no --known_sites_index and engine options..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  "${reference_args[@]}" \
  --known_sites "$test_data_dir/known_sites1.vcf.gz" \
  --interval_padding 10 \
  --read_filter MappingQualityReadFilter \
  --disable_tool_default_read_filters \
  --output "$meta_temp_dir/recal_no_index.table" 2> "$meta_temp_dir/no_index.log"

check_file_contains "$meta_temp_dir/recal_no_index.table" "#:GATKReport" "recalibration table GATKReport header (no known-sites index)"
check_file_contains "$meta_temp_dir/no_index.log" "no index was provided" "warning for the missing index"
check_file_not_exists "$test_data_dir/known_sites1.vcf.gz.tbi" "index next to the known-sites input"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Multiple --intervals files (BED and bgzipped VCF) with --interval_set_rule INTERSECTION"

# Reads start at 1, 101, ..., 901 (100bp each). The area is the
# intersection of seq1:401-800 (BED) and seq1:650 (VCF), so only the read at
# 601 is used (100 bases, no known site in it). With UNION, the reads at 401,
# 501, 601 and 701 would be used.
printf 'seq1\t400\t800\n' > "$test_data_dir/intervals_a.bed"
create_test_known_sites_vcf --no-index "$test_data_dir/intervals_site.vcf" seq1 1000 650
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  "${reference_args[@]}" \
  --known_sites "$known_sites1" \
  --output "$meta_temp_dir/recal_intervals.table" \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule INTERSECTION

log "Validating TEST 3 outputs..."
observations=$(rg_observations "$meta_temp_dir/recal_intervals.table" rg1)
if [[ "$observations" == "100" ]]; then
  log "✓ recalibration table only has the bases of the read in the selected territory ($observations observations)"
else
  log_error "✗ recalibration table has '$observations' observations, expected 100"
  exit 1
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Argument errors ---
log "Starting TEST 4: Argument errors"

log "Executing $meta_name with two known-sites files and one index (should fail)..."
check_fails_with "Error: --known_sites_index must be given once for each file in --known_sites" "$meta_temp_dir/index_count.log" \
  "$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --read_index "$test_data_dir/reads.sorted.bai" \
  "${reference_args[@]}" \
  --known_sites "$test_data_dir/known_sites1.vcf.gz" \
  --known_sites "$test_data_dir/known_sites2.vcf.gz" \
  --known_sites_index "$test_data_dir/indexes/known_sites1.vcf.gz.tbi" \
  --output "$meta_temp_dir/recal_bad.table"
check_file_not_exists "$meta_temp_dir/recal_bad.table" "recalibration table (should not be written)"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
