#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

# Two per-sample GVCFs (plain VCF) from a direct HaplotypeCaller run
for sample in sample1 sample2; do
  create_test_sam_reads "$test_data_dir/$sample.sam" seq1 2000 "rg_$sample" "$sample"
  sort_and_index_bam "$test_data_dir/$sample.sam" "$test_data_dir/$sample.sorted.bam"
  check_file_exists "$test_data_dir/$sample.sorted.bai" "BAM index file for $sample"

  log "Running gatk HaplotypeCaller on $sample..."
  gatk HaplotypeCaller \
    --input "$test_data_dir/$sample.sorted.bam" \
    --reference "$test_data_dir/reference.fasta" \
    --emit-ref-confidence GVCF \
    --output "$test_data_dir/$sample.g.vcf" \
    --verbosity ERROR
  check_file_not_empty "$test_data_dir/$sample.g.vcf" "$sample GVCF"
done

# An intervals file that covers the reference contig
echo "seq1:1-2000" > "$test_data_dir/intervals.list"

# A sample name map for the plain GVCFs
printf 'sample1\t%s\nsample2\t%s\n' \
  "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/sample2.g.vcf" \
  > "$test_data_dir/sample_map.tsv"

# Bgzipped GVCFs (and a header VCF) with their indexes in a different
# directory
create_test_vcf_gz_with_separate_index "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/separate/sample1.g.vcf.gz" "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sample2.g.vcf" \
  "$test_data_dir/separate/sample2.g.vcf.gz" "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/separate/header.vcf.gz" "$test_data_dir/separate_index/header.vcf.gz.tbi"

# Bgzipped GVCFs (and a header VCF) without indexes
mkdir -p "$test_data_dir/no_index"
bgzip -c "$test_data_dir/sample1.g.vcf" > "$test_data_dir/no_index/sample1.g.vcf.gz"
bgzip -c "$test_data_dir/sample2.g.vcf" > "$test_data_dir/no_index/sample2.g.vcf.gz"
bgzip -c "$test_data_dir/sample1.g.vcf" > "$test_data_dir/no_index/header.vcf.gz"

# Bgzipped GVCFs with their indexes next to them, and a sample name map for
# them (--bypass_feature_reader and --avoid_nio need block-compressed and
# indexed GVCFs)
mkdir -p "$test_data_dir/bgzip"
for sample in sample1 sample2; do
  bgzip -c "$test_data_dir/$sample.g.vcf" > "$test_data_dir/bgzip/$sample.g.vcf.gz"
  tabix -p vcf "$test_data_dir/bgzip/$sample.g.vcf.gz"
  check_file_exists "$test_data_dir/bgzip/$sample.g.vcf.gz.tbi" "$sample block-compressed GVCF index"
done
printf 'sample1\t%s\nsample2\t%s\n' \
  "$test_data_dir/bgzip/sample1.g.vcf.gz" \
  "$test_data_dir/bgzip/sample2.g.vcf.gz" \
  > "$test_data_dir/sample_map_gz.tsv"

# Interval files: include seq1:100-200 (BED is 0-based) and seq1:1001-1020
# (the 20bp REF allele of the VCF record), exclude seq1:150-160
mkdir -p "$test_data_dir/intervals"
printf 'seq1\t99\t200\n' > "$test_data_dir/intervals/include.bed"
printf 'seq1\t149\t160\n' > "$test_data_dir/intervals/exclude.bed"
{
  echo "##fileformat=VCFv4.2"
  echo "##contig=<ID=seq1,length=2000>"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'seq1\t1001\t.\tATCGATCGATCGATCGATCG\tA\t50\tPASS\t.\n'
} > "$test_data_dir/intervals/include.vcf"
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals/include.vcf" \
  "$test_data_dir/intervals/include.vcf.gz" "$test_data_dir/intervals_index/include.vcf.gz.tbi"

# Print the GenomicsDB array directories (one for each imported interval,
# named <contig>$<start>$<end>) in a workspace
list_arrays() {
  find "$1" -mindepth 1 -maxdepth 1 -type d -name 'seq1$*' -exec basename {} \; | LC_ALL=C sort | tr '\n' ' '
}

# Check that a workspace has exactly the expected arrays
#
# Usage: check_arrays "/path/to/workspace" "expected arrays " "description"
check_arrays() {
  local arrays
  arrays=$(list_arrays "$1")
  if [[ "$arrays" == "$2" ]]; then
    log "✓ $3: $arrays"
  else
    log_error "✗ Workspace arrays are '$arrays', expected '$2' ($3)"
    exit 1
  fi
}

# --- Test Case 1: Basic run ---
log "Starting TEST 1: Basic GenomicsDBImport with --variant"

"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace"

log "Validating TEST 1 outputs..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace" "GenomicsDB workspace directory"
check_file_exists "$meta_temp_dir/genomicsdb_workspace/callset.json" "GenomicsDB callset.json"
check_file_exists "$meta_temp_dir/genomicsdb_workspace/vidmap.json" "GenomicsDB vidmap.json"
check_file_exists "$meta_temp_dir/genomicsdb_workspace/vcfheader.vcf" "GenomicsDB vcfheader.vcf"
check_file_contains "$meta_temp_dir/genomicsdb_workspace/callset.json" "sample1" "GenomicsDB callset.json sample1"
check_file_contains "$meta_temp_dir/genomicsdb_workspace/callset.json" "sample2" "GenomicsDB callset.json sample2"

# Confirm the per-contig array subdirectory was created inside the workspace
contig_subdir_count=$(find "$meta_temp_dir/genomicsdb_workspace" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
if [[ "$contig_subdir_count" -lt 1 ]]; then
  log_error "Expected at least one per-contig array subdirectory in the GenomicsDB workspace"
  exit 1
fi

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Input index handling, --sample_name, --sample_name_map and --header"

# Run 1: bgzipped --variant files with separate indexes (no --sample_name, so
# GATK gets --variant)
log "Executing $meta_name with .vcf.gz --variant files and separate .tbi index files..."
"$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --variant_index "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_separate" 2> "$meta_temp_dir/test2_separate.log"

log "Validating TEST 2 outputs (separate --variant indexes)..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_separate/callset.json" "GenomicsDB callset.json (separate indexes)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_separate/callset.json" "sample1" "GenomicsDB callset.json sample1 (separate indexes)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_separate/callset.json" "sample2" "GenomicsDB callset.json sample2 (separate indexes)"
check_file_not_contains "$meta_temp_dir/test2_separate.log" "no index was provided" "run log (no indexing when indexes are given)"
check_file_not_exists "$test_data_dir/separate/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/separate/sample2.g.vcf.gz.tbi" "index next to the sample2 input"

# Run 2: unindexed bgzipped --variant files renamed with --sample_name, and an
# unindexed bgzipped --header
log "Executing $meta_name with unindexed .vcf.gz --variant files, --sample_name and an unindexed .vcf.gz --header..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/sample1.g.vcf.gz" \
  --variant "$test_data_dir/no_index/sample2.g.vcf.gz" \
  --sample_name renamed1 \
  --sample_name renamed2 \
  --header "$test_data_dir/no_index/header.vcf.gz" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_no_index" 2> "$meta_temp_dir/test2_no_index.log"

log "Validating TEST 2 outputs (unindexed inputs, --sample_name)..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_no_index/vcfheader.vcf" "GenomicsDB vcfheader.vcf (no header index)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_no_index/callset.json" "renamed1" "GenomicsDB callset.json renamed1"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_no_index/callset.json" "renamed2" "GenomicsDB callset.json renamed2"
check_file_not_contains "$meta_temp_dir/genomicsdb_workspace_no_index/callset.json" "\"sample1\"" "GenomicsDB callset.json (original sample1 name)"
check_file_contains "$meta_temp_dir/test2_no_index.log" "no index was provided for '.*no_index/sample1.g.vcf.gz'" "run log (indexing warning for --variant)"
check_file_contains "$meta_temp_dir/test2_no_index.log" "no index was provided for '.*no_index/header.vcf.gz'" "run log (indexing warning for --header)"
check_file_not_exists "$test_data_dir/no_index/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/no_index/header.vcf.gz.tbi" "index next to the header input"

# Run 3: --sample_name_map with the import options, and a bgzipped --header
# with a separate index
log "Executing $meta_name with --sample_name_map, import options and a .vcf.gz --header with --header_index..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_map" \
  --header "$test_data_dir/separate/header.vcf.gz" \
  --header_index "$test_data_dir/separate_index/header.vcf.gz.tbi" \
  --batch_size 1 \
  --consolidate \
  --validate_sample_name_map \
  --genomicsdb_vcf_buffer_size 32768 \
  --max_num_intervals_to_import_in_parallel 1 \
  --merge_contigs_into_num_partitions 1 > "$meta_temp_dir/test2_map.log" 2>&1

log "Validating TEST 2 outputs (--sample_name_map, --header_index)..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace_map" "GenomicsDB workspace directory (sample_name_map)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_map/callset.json" "GenomicsDB callset.json (sample_name_map)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_map/vidmap.json" "GenomicsDB vidmap.json (sample_name_map)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_map/vcfheader.vcf" "GenomicsDB vcfheader.vcf (sample_name_map)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_map/callset.json" "sample1" "GenomicsDB callset.json sample1 (sample_name_map)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_map/callset.json" "sample2" "GenomicsDB callset.json sample2 (sample_name_map)"
check_file_not_contains "$meta_temp_dir/test2_map.log" "no index was provided" "run log (no indexing when an index is given)"
check_file_not_exists "$test_data_dir/separate/header.vcf.gz.tbi" "index next to the header input"

for expected in \
  "--batch-size 1" \
  "--consolidate" \
  "--validate-sample-name-map" \
  "--genomicsdb-vcf-buffer-size 32768" \
  "--max-num-intervals-to-import-in-parallel 1" \
  "--merge-contigs-into-num-partitions 1"; do
  check_file_contains "$meta_temp_dir/test2_map.log" "GenomicsDBImport .* $expected" "GATK command line ($expected)"
done

# Run 4: --bypass_feature_reader and --avoid_nio. --avoid_nio cannot be used
# together with --variant, only with --sample_name_map, and GATK also requires
# --header whenever --avoid_nio is used.
log "Executing $meta_name with --bypass_feature_reader and --avoid_nio..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map_gz.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_bypass" \
  --bypass_feature_reader \
  --avoid_nio \
  --header "$test_data_dir/bgzip/sample1.g.vcf.gz"

log "Validating TEST 2 outputs (--bypass_feature_reader, --avoid_nio)..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace_bypass" "GenomicsDB workspace directory (bypass_feature_reader)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_bypass/callset.json" "GenomicsDB callset.json (bypass_feature_reader)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_bypass/callset.json" "sample1" "GenomicsDB callset.json sample1 (bypass_feature_reader)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_bypass/callset.json" "sample2" "GenomicsDB callset.json sample2 (bypass_feature_reader)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Interval options"

# Run 1: only --interval_strings, with the engine options. --interval_padding
# 10 extends seq1:100-200 to seq1:90-210.
log "Executing $meta_name with only --interval_strings and engine options..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --interval_strings "seq1:100-200" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_strings" \
  --interval_padding 10 \
  --lenient \
  --disable_sequence_dictionary_validation > "$meta_temp_dir/test3_strings.log" 2>&1

log "Validating TEST 3 outputs (--interval_strings, engine options)..."
check_file_contains "$meta_temp_dir/genomicsdb_workspace_strings/callset.json" "sample1" "GenomicsDB callset.json sample1 (engine options)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_strings/callset.json" "sample2" "GenomicsDB callset.json sample2 (engine options)"
check_arrays "$meta_temp_dir/genomicsdb_workspace_strings" 'seq1$90$210 ' "Workspace has one array for the padded interval string"
for expected in "--lenient" "--disable-sequence-dictionary-validation"; do
  check_file_contains "$meta_temp_dir/test3_strings.log" "GenomicsDBImport .* $expected" "GATK command line ($expected)"
done

# Run 2: two --intervals files and --exclude_intervals
log "Executing $meta_name with two --intervals files and --exclude_intervals..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --intervals "$test_data_dir/intervals/include.bed" \
  --intervals "$test_data_dir/intervals/include.vcf.gz" \
  --intervals_index "$test_data_dir/intervals_index/include.vcf.gz.tbi" \
  --exclude_intervals "$test_data_dir/intervals/exclude.bed" \
  --interval_set_rule UNION \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_intervals" 2> "$meta_temp_dir/test3_intervals.log"

log "Validating TEST 3 outputs (interval files)..."
check_file_contains "$meta_temp_dir/genomicsdb_workspace_intervals/callset.json" "sample1" "GenomicsDB callset.json sample1 (intervals)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_intervals/callset.json" "sample2" "GenomicsDB callset.json sample2 (intervals)"
check_file_not_contains "$meta_temp_dir/test3_intervals.log" "no index was provided for '$test_data_dir/intervals/include.vcf.gz'" "run log (no indexing of the interval VCF when --intervals_index is given)"
check_file_not_exists "$test_data_dir/intervals/include.vcf.gz.tbi" "index next to the interval VCF"
check_arrays "$meta_temp_dir/genomicsdb_workspace_intervals" 'seq1$100$149 seq1$1001$1020 seq1$161$200 ' "Workspace has one array for each selected interval"

# Run 3: two --intervals files and --merge_input_intervals
log "Executing $meta_name with two --intervals files and --merge_input_intervals..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --intervals "$test_data_dir/intervals/include.bed" \
  --intervals "$test_data_dir/intervals/include.vcf.gz" \
  --intervals_index "$test_data_dir/intervals_index/include.vcf.gz.tbi" \
  --merge_input_intervals \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_merged"

log "Validating TEST 3 outputs (--merge_input_intervals)..."
check_arrays "$meta_temp_dir/genomicsdb_workspace_merged" 'seq1$100$1020 ' "Workspace has one merged array for the contig"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Argument errors ---
log "Starting TEST 4: Argument errors"

log "Executing $meta_name with 2 --variant files and 1 --variant_index file..."
check_fails_with "Error: --variant_index must be given once for each file in --variant" "$meta_temp_dir/test4_mismatch.log" \
  "$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_mismatch"

log "Executing $meta_name with --variant_index and no --variant..."
check_fails_with "Error: --variant_index requires --variant" "$meta_temp_dir/test4_variant.log" \
  "$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_variant_index"

log "Executing $meta_name with --header_index and no --header..."
check_fails_with "Error: --header_index requires --header" "$meta_temp_dir/test4_header.log" \
  "$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --header_index "$test_data_dir/separate_index/header.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_header_index"

log "Executing $meta_name with --header and --variant..."
check_fails_with "Error: --header with --variant requires --sample_name" "$meta_temp_dir/test4_header_variant.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --header "$test_data_dir/sample1.g.vcf" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_header_variant"

log "Executing $meta_name with fewer --sample_name values than --variant files..."
check_fails_with "Error: --sample_name must be given once for each file in --variant" "$meta_temp_dir/test4_sample_name_count.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --sample_name renamed1 \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_sample_name_count"

log "Executing $meta_name with --sample_name and --sample_name_map..."
check_fails_with "Error: --sample_name cannot be used together with --sample_name_map" "$meta_temp_dir/test4_sample_name_map.log" \
  "$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --sample_name renamed1 \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_sample_name_map"

log "Executing $meta_name without --intervals or --interval_strings..."
check_fails_with "Error: --intervals or --interval_strings is required." "$meta_temp_dir/test4_no_intervals.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_no_intervals"
check_dir_not_exists "$meta_temp_dir/genomicsdb_workspace_no_intervals" "GenomicsDB workspace (no run without --intervals)"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
