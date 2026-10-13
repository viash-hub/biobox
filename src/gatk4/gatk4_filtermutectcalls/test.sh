#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Set up the test environment and the test data directory
gatk4_test_init

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

# Regex prefix for the chr1:250 call up to its FILTER column (column 7)
call_250_filter=$'^chr1\t250\t([^\t]*\t){4}'

# --- Test Case 1: Required inputs only ---
log "Starting TEST 1: Required inputs only"

log "Executing $meta_name with an unindexed --variant and --stats..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/filtered.vcf" \
  --output_index "$meta_temp_dir/filtered.idx" \
  --filtering_stats "$meta_temp_dir/filtering_stats.tsv" 2>&1 | tee "$meta_temp_dir/test1.log"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/filtered.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/filtered.vcf" "^##fileformat=VCF" "output VCF file header"
check_vcf_has_record "$meta_temp_dir/filtered.vcf" "${call_250_filter}PASS"$'\t' "output VCF (chr1:250 call has FILTER PASS)"
check_file_exists "$meta_temp_dir/filtering_stats.tsv" "output filtering stats file"
check_file_not_empty "$meta_temp_dir/filtering_stats.tsv" "output filtering stats file"
check_file_not_exists "$meta_temp_dir/filtered.vcf.filteringStats.tsv" "default filtering stats file next to the output VCF"
check_file_not_contains "$meta_temp_dir/test1.log" "no index was provided" "run log (plain --variant used without an index)"
check_file_not_exists "$test_data_dir/no_index/unfiltered.vcf.idx" "index next to the plain --variant"

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
check_vcf_has_record "$meta_temp_dir/filtered_all.vcf.gz" "${call_250_filter}"$'([^.\t][^\t]*|\\.[^\t]+)\t' "output VCF (chr1:250 call has a FILTER value)"
check_index_moved "$meta_temp_dir/filtered_all.tbi" "$meta_temp_dir/filtered_all.vcf.gz.tbi" "output index file"
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
  --output_index "$meta_temp_dir/filtered_strand.idx" \
  --filtering_stats "$meta_temp_dir/filtering_stats_strand.tsv"

log "Validating TEST 3 outputs..."
check_vcf_has_record "$meta_temp_dir/filtered_strand.vcf" "${call_250_filter}"$'[^\t]*strict_strand' "output VCF (chr1:250 call has the strict_strand filter)"
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
check_vcf_record_count "$meta_temp_dir/filtered_intervals.vcf" 0 "output VCF (no calls in --intervals)"
check_file_not_exists "$meta_temp_dir/filtered_intervals.idx" "output index file"
check_file_not_exists "$meta_temp_dir/filtered_intervals.vcf.idx" "index next to the output VCF"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Missing --output_index ---
log "Starting TEST 5: A missing --output_index fails"

check_fails_with "Error: --output_index is required unless" "$meta_temp_dir/test5.log" \
  "$meta_executable" \
  --variant "$test_data_dir/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_output_index.vcf" \
  --filtering_stats "$meta_temp_dir/no_output_index.tsv"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Interval padding ---
log "Starting TEST 6: --interval_padding and --interval_exclusion_padding"

# Positions 241-245, which end 5bp before the call at position 250
printf 'chr1\t240\t245\n' > "$test_data_dir/near_250.bed"

log "Executing $meta_name with an --intervals file next to the call, without padding..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --intervals "$test_data_dir/near_250.bed" \
  --output "$meta_temp_dir/filtered_no_padding.vcf" \
  --output_index "$meta_temp_dir/filtered_no_padding.idx" \
  --filtering_stats "$meta_temp_dir/filtering_stats_no_padding.tsv"

check_vcf_lacks_record "$meta_temp_dir/filtered_no_padding.vcf" $'^chr1\t250\t' "output VCF (chr1:250 call outside the interval)"

log "Executing $meta_name with the same --intervals file and --interval_padding 5..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --intervals "$test_data_dir/near_250.bed" \
  --interval_padding 5 \
  --output "$meta_temp_dir/filtered_padding.vcf" \
  --output_index "$meta_temp_dir/filtered_padding.idx" \
  --filtering_stats "$meta_temp_dir/filtering_stats_padding.tsv"

check_vcf_has_record "$meta_temp_dir/filtered_padding.vcf" "${call_250_filter}PASS"$'\t' "output VCF (chr1:250 call in the padded interval, FILTER PASS)"

# Only --exclude_intervals is given, so the plain --variant must also be
# indexed (GATK traverses by intervals)
log "Executing $meta_name with the same file as --exclude_intervals and --interval_exclusion_padding 5..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/unfiltered.vcf" \
  --stats "$test_data_dir/unfiltered.vcf.stats" \
  "${reference_args[@]}" \
  --exclude_intervals "$test_data_dir/near_250.bed" \
  --interval_exclusion_padding 5 \
  --output "$meta_temp_dir/filtered_exclusion_padding.vcf" \
  --output_index "$meta_temp_dir/filtered_exclusion_padding.idx" \
  --filtering_stats "$meta_temp_dir/filtering_stats_exclusion_padding.tsv"

check_vcf_lacks_record "$meta_temp_dir/filtered_exclusion_padding.vcf" $'^chr1\t250\t' "output VCF (chr1:250 call in the padded excluded interval)"
check_file_not_exists "$test_data_dir/no_index/unfiltered.vcf.idx" "index next to the plain --variant"

log "✅ TEST 6 completed successfully"

print_test_summary "All tests"
