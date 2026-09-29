#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
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

# --- Build shared fixtures: a real Mutect2 run on a tumor/normal pair ---
create_test_mutect2_outputs "$test_data_dir"
check_file_exists "$test_data_dir/unfiltered.vcf" "unfiltered Mutect2 VCF"
check_file_exists "$test_data_dir/unfiltered.vcf.stats" "Mutect2 stats file"

# A copy of the VCF without an index next to it
mkdir -p "$test_data_dir/no_index"
cp "$test_data_dir/unfiltered.vcf" "$test_data_dir/no_index/unfiltered.vcf"

# A bgzipped copy of the VCF with a tabix index
mkdir -p "$test_data_dir/indexed"
bgzip -c "$test_data_dir/unfiltered.vcf" > "$test_data_dir/indexed/unfiltered.vcf.gz"
tabix -p vcf "$test_data_dir/indexed/unfiltered.vcf.gz"

reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

# Print the FILTER value of the chr1:250 call in a VCF
filter_at_250() {
  local vcf_path="$1"
  if [[ "$vcf_path" == *.gz ]]; then
    zcat "$vcf_path"
  else
    cat "$vcf_path"
  fi | awk -F '\t' '!/^#/ && $1 == "chr1" && $2 == 250 { print $7 }'
}

# --- Test Case 1: Required inputs only ---
log "Starting TEST 1: Required inputs only"

log "Executing $meta_name with an unindexed --variant and --stats..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/filtered.vcf" \
  --filtering_stats "$meta_temp_dir/filtering_stats.tsv" 2>&1 | tee "$meta_temp_dir/test1.log"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/filtered.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/filtered.vcf" "^##fileformat=VCF" "output VCF file header"
filter_value=$(filter_at_250 "$meta_temp_dir/filtered.vcf")
if [[ "$filter_value" == "PASS" ]]; then
  log "✓ The chr1:250 call has FILTER PASS"
else
  log_error "✗ Expected FILTER PASS for the chr1:250 call, found '$filter_value'"
  exit 1
fi
check_file_exists "$meta_temp_dir/filtering_stats.tsv" "output filtering stats file"
check_file_not_empty "$meta_temp_dir/filtering_stats.tsv" "output filtering stats file"
check_file_not_exists "$meta_temp_dir/filtered.vcf.filteringStats.tsv" "default filtering stats file next to the output VCF"
check_file_contains "$meta_temp_dir/test1.log" "Warning: no index was provided for '.*no_index/unfiltered.vcf'" "run log (copy warning for the unindexed --variant)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: All optional filtering inputs and an output index ---
log "Starting TEST 2: All optional filtering inputs and an output index"

log "Writing contamination and tumor segmentation tables..."
printf 'sample\tcontamination\terror\ntumor\t0.02\t0.01\n' > "$test_data_dir/contamination.table"
printf 'contig\tstart\tend\tminor_allele_fraction\nchr1\t1\t500\t0.5\n' > "$test_data_dir/segmentation.table"

log "Learning a read orientation model from the F1R2 counts..."
gatk LearnReadOrientationModel \
  --input "$test_data_dir/f1r2.tar.gz" \
  --output "$test_data_dir/read_orientation_model.tar.gz" \
  --verbosity ERROR
check_file_exists "$test_data_dir/read_orientation_model.tar.gz" "read orientation model"

log "Executing $meta_name with an indexed .vcf.gz --variant and all optional filtering inputs..."
"$meta_executable" \
  --variant "$test_data_dir/indexed/unfiltered.vcf.gz" \
  --variant_index "$test_data_dir/indexed/unfiltered.vcf.gz.tbi" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --contamination_table "$test_data_dir/contamination.table" \
  --tumor_segmentation "$test_data_dir/segmentation.table" \
  --orientation_bias_artifact_priors "$test_data_dir/read_orientation_model.tar.gz" \
  --output "$meta_temp_dir/filtered_all.vcf.gz" \
  --output_index "$meta_temp_dir/filtered_all.tbi" \
  --filtering_stats "$meta_temp_dir/filtering_stats_all.tsv" 2>&1 | tee "$meta_temp_dir/test2.log"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/filtered_all.vcf.gz" "output VCF file"
filter_value=$(filter_at_250 "$meta_temp_dir/filtered_all.vcf.gz")
if [[ -n "$filter_value" && "$filter_value" != "." ]]; then
  log "✓ The chr1:250 call has a FILTER value: $filter_value"
else
  log_error "✗ The chr1:250 call has no FILTER value: '$filter_value'"
  exit 1
fi
check_file_exists "$meta_temp_dir/filtered_all.tbi" "output index file"
check_file_not_exists "$meta_temp_dir/filtered_all.vcf.gz.tbi" "index next to the output VCF (moved to --output_index)"
check_file_exists "$meta_temp_dir/filtering_stats_all.tsv" "output filtering stats file"
check_file_not_contains "$meta_temp_dir/test2.log" "Warning: no index was provided" "run log (no copy warning when the index is given)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Threshold and hard filter options ---
log "Starting TEST 3: Threshold and hard filter options"

# The alternate reads of the fixture are all on the forward strand, so a
# minimum of 1 alternate read per strand filters the call
log "Executing $meta_name with --threshold_strategy CONSTANT and --min_reads_per_strand 1..."
"$meta_executable" \
  --variant "$test_data_dir/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --threshold_strategy CONSTANT \
  --initial_threshold 0.1 \
  --min_reads_per_strand 1 \
  --output "$meta_temp_dir/filtered_strand.vcf" \
  --filtering_stats "$meta_temp_dir/filtering_stats_strand.tsv"

log "Validating TEST 3 outputs..."
filter_value=$(filter_at_250 "$meta_temp_dir/filtered_strand.vcf")
if [[ "$filter_value" == *strict_strand* ]]; then
  log "✓ The chr1:250 call has the strict_strand filter: $filter_value"
else
  log_error "✗ Expected the strict_strand FILTER for the chr1:250 call, found '$filter_value'"
  exit 1
fi
check_file_contains "$meta_temp_dir/filtered_strand.vcf" "threshold-strategy CONSTANT" "output VCF command line (--threshold_strategy)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Intervals and output index without index creation ---
log "Starting TEST 4: Intervals, and --output_index without index creation"

printf 'chr1\t0\t100\n' > "$test_data_dir/no_calls.bed"

# The --output_index path is set with --create_output_variant_index false, the
# same as the Nextflow runner does, so it must be ignored
log "Executing $meta_name with an --intervals BED file that has no calls..."
"$meta_executable" \
  --variant "$test_data_dir/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --intervals "$test_data_dir/no_calls.bed" \
  --output "$meta_temp_dir/filtered_intervals.vcf" \
  --output_index "$meta_temp_dir/filtered_intervals.idx" \
  --create_output_variant_index false \
  --filtering_stats "$meta_temp_dir/filtering_stats_intervals.tsv"

log "Validating TEST 4 outputs..."
record_count=$(grep -vc '^#' "$meta_temp_dir/filtered_intervals.vcf" || true)
if [[ "$record_count" -eq 0 ]]; then
  log "✓ No calls outside --intervals"
else
  log_error "✗ Expected no calls in --intervals, found $record_count"
  exit 1
fi
check_file_not_exists "$meta_temp_dir/filtered_intervals.idx" "output index file"
check_file_not_exists "$meta_temp_dir/filtered_intervals.vcf.idx" "index next to the output VCF"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
