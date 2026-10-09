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

# --- Build shared reference fixture ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

# --- Build two single-sample BAM fixtures with distinct read groups ---
create_test_sam_reads "$test_data_dir/sample1.sam" seq1 2000 rg1 sample1
check_file_exists "$test_data_dir/sample1.sam" "hand-crafted SAM reads for sample1"

create_test_sam_reads "$test_data_dir/sample2.sam" seq1 2000 rg2 sample2
check_file_exists "$test_data_dir/sample2.sam" "hand-crafted SAM reads for sample2"

log "Sorting and indexing sample1 reads..."
sort_and_index_bam "$test_data_dir/sample1.sam" "$test_data_dir/sample1.sorted.bam"
check_file_exists "$test_data_dir/sample1.sorted.bam" "sorted BAM file for sample1"
check_file_exists "$test_data_dir/sample1.sorted.bai" "BAM index file for sample1"

log "Sorting and indexing sample2 reads..."
sort_and_index_bam "$test_data_dir/sample2.sam" "$test_data_dir/sample2.sorted.bam"
check_file_exists "$test_data_dir/sample2.sorted.bam" "sorted BAM file for sample2"
check_file_exists "$test_data_dir/sample2.sorted.bai" "BAM index file for sample2"

# --- Generate per-sample GVCFs directly via gatk HaplotypeCaller ---
log "Running gatk HaplotypeCaller on sample1..."
gatk HaplotypeCaller \
  --input "$test_data_dir/sample1.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --output "$test_data_dir/sample1.g.vcf" \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample1.g.vcf" "sample1 GVCF"
check_file_not_empty "$test_data_dir/sample1.g.vcf" "sample1 GVCF"

log "Running gatk HaplotypeCaller on sample2..."
gatk HaplotypeCaller \
  --input "$test_data_dir/sample2.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --output "$test_data_dir/sample2.g.vcf" \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample2.g.vcf" "sample2 GVCF"
check_file_not_empty "$test_data_dir/sample2.g.vcf" "sample2 GVCF"

# --- Build intervals file covering the reference contig ---
log "Writing intervals file..."
cat > "$test_data_dir/intervals.list" << 'EOF'
seq1:1-2000
EOF
check_file_exists "$test_data_dir/intervals.list" "intervals file"

# --- Test Case 1: Basic import using repeated --variant flags ---
log "Starting TEST 1: Basic GenomicsDBImport with --variant"

log "Executing $meta_name with basic parameters..."
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

# --- Test Case 2: Import using --sample_name_map instead of --variant ---
log "Starting TEST 2: GenomicsDBImport with --sample_name_map"

log "Writing sample name map..."
printf 'sample1\t%s\nsample2\t%s\n' \
  "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/sample2.g.vcf" \
  > "$test_data_dir/sample_map.tsv"
check_file_exists "$test_data_dir/sample_map.tsv" "sample name map"

log "Executing $meta_name with --sample_name_map and extra options..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_map" \
  --batch_size 1 \
  --consolidate \
  --validate_sample_name_map \
  --genomicsdb_vcf_buffer_size 32768 \
  --max_num_intervals_to_import_in_parallel 1 \
  --merge_contigs_into_num_partitions 1

log "Validating TEST 2 outputs..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace_map" "GenomicsDB workspace directory (sample_name_map)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_map/callset.json" "GenomicsDB callset.json (sample_name_map)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_map/vidmap.json" "GenomicsDB vidmap.json (sample_name_map)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_map/vcfheader.vcf" "GenomicsDB vcfheader.vcf (sample_name_map)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_map/callset.json" "sample1" "GenomicsDB callset.json sample1 (sample_name_map)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_map/callset.json" "sample2" "GenomicsDB callset.json sample2 (sample_name_map)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Import using --header to supply a fixed VCF header ---
log "Starting TEST 3: GenomicsDBImport with --header"

log "Executing $meta_name with --header..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_header" \
  --header "$test_data_dir/sample1.g.vcf"

log "Validating TEST 3 outputs..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace_header" "GenomicsDB workspace directory (header)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_header/callset.json" "GenomicsDB callset.json (header)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_header/vcfheader.vcf" "GenomicsDB vcfheader.vcf (header)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_header/callset.json" "sample1" "GenomicsDB callset.json sample1 (header)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_header/callset.json" "sample2" "GenomicsDB callset.json sample2 (header)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Import using --bypass_feature_reader/--avoid_nio ---
# These options require normalized, block-compressed and indexed GVCFs, so
# generate block-gzipped GVCFs directly via HaplotypeCaller (writing to a
# ".g.vcf.gz" output makes GATK/htsjdk emit a block-compressed VCF with a
# .tbi index, without depending on external bgzip/tabix binaries).
log "Starting TEST 4: GenomicsDBImport with --bypass_feature_reader/--avoid_nio"

log "Running gatk HaplotypeCaller on sample1 (block-compressed output)..."
gatk HaplotypeCaller \
  --input "$test_data_dir/sample1.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --output "$test_data_dir/sample1.g.vcf.gz" \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample1.g.vcf.gz" "sample1 block-compressed GVCF"
check_file_exists "$test_data_dir/sample1.g.vcf.gz.tbi" "sample1 block-compressed GVCF index"

log "Running gatk HaplotypeCaller on sample2 (block-compressed output)..."
gatk HaplotypeCaller \
  --input "$test_data_dir/sample2.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --output "$test_data_dir/sample2.g.vcf.gz" \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample2.g.vcf.gz" "sample2 block-compressed GVCF"
check_file_exists "$test_data_dir/sample2.g.vcf.gz.tbi" "sample2 block-compressed GVCF index"

log "Writing sample name map for block-compressed GVCFs..."
printf 'sample1\t%s\nsample2\t%s\n' \
  "$test_data_dir/sample1.g.vcf.gz" \
  "$test_data_dir/sample2.g.vcf.gz" \
  > "$test_data_dir/sample_map_gz.tsv"
check_file_exists "$test_data_dir/sample_map_gz.tsv" "sample name map (block-compressed)"

# Note: --avoid_nio cannot be used together with --variant, only with
# --sample_name_map, and GATK additionally requires --header to be set
# whenever --avoid_nio is used
log "Executing $meta_name with --bypass_feature_reader and --avoid_nio..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map_gz.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_bypass" \
  --bypass_feature_reader \
  --avoid_nio \
  --header "$test_data_dir/sample1.g.vcf.gz"

log "Validating TEST 4 outputs..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace_bypass" "GenomicsDB workspace directory (bypass_feature_reader)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_bypass/callset.json" "GenomicsDB callset.json (bypass_feature_reader)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_bypass/callset.json" "sample1" "GenomicsDB callset.json sample1 (bypass_feature_reader)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_bypass/callset.json" "sample2" "GenomicsDB callset.json sample2 (bypass_feature_reader)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Shared GATK engine options ---
# GenomicsDBImport's output is a workspace directory, not a VCF/BAM file, so
# --sites_only_vcf_output/--create_output_variant_index/--create_output_bam_index
# are semantically inert here, but GATK accepts them without error, so we
# just confirm the tool still runs successfully with all the flags set
log "Starting TEST 5: Shared GATK engine options"

log "Executing $meta_name with --interval_padding and other shared engine options..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_engine" \
  --interval_padding 10 \
  --sites_only_vcf_output \
  --create_output_variant_index false \
  --create_output_bam_index false \
  --disable_tool_default_read_filters \
  --disable_sequence_dictionary_validation

log "Validating TEST 5 outputs..."
check_dir_exists "$meta_temp_dir/genomicsdb_workspace_engine" "GenomicsDB workspace directory (engine options)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_engine/callset.json" "GenomicsDB callset.json (engine options)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_engine/vidmap.json" "GenomicsDB vidmap.json (engine options)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_engine/vcfheader.vcf" "GenomicsDB vcfheader.vcf (engine options)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_engine/callset.json" "sample1" "GenomicsDB callset.json sample1 (engine options)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_engine/callset.json" "sample2" "GenomicsDB callset.json sample2 (engine options)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Index files in a different directory ---
# In GATK 4.7.0.0, --header together with --variant fails with
# "java.lang.ArithmeticException: / by zero", so --header is tested with
# --sample_name_map
log "Starting TEST 6: --variant_index and --header_index in a different directory"

mkdir -p "$test_data_dir/separate"

create_test_vcf_gz_with_separate_index "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/separate/sample1.g.vcf.gz" "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sample2.g.vcf" \
  "$test_data_dir/separate/sample2.g.vcf.gz" "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/separate/header.vcf.gz" "$test_data_dir/separate_index/header.vcf.gz.tbi"

log "Executing $meta_name with .vcf.gz --variant files and separate .tbi index files..."
"$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --variant_index "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_separate" 2> "$meta_temp_dir/test6.log"

log "Validating TEST 6 --variant outputs..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_separate/callset.json" "GenomicsDB callset.json (separate indexes)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_separate/callset.json" "sample1" "GenomicsDB callset.json sample1 (separate indexes)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_separate/callset.json" "sample2" "GenomicsDB callset.json sample2 (separate indexes)"
check_file_not_contains "$meta_temp_dir/test6.log" "no index was provided" "run log (no indexing when indexes are given)"
check_file_not_exists "$test_data_dir/separate/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/separate/sample2.g.vcf.gz.tbi" "index next to the sample2 input"

log "Executing $meta_name with --sample_name_map, a .vcf.gz --header and a separate .tbi index file..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_separate_header" \
  --header "$test_data_dir/separate/header.vcf.gz" \
  --header_index "$test_data_dir/separate_index/header.vcf.gz.tbi" 2> "$meta_temp_dir/test6_header.log"

log "Validating TEST 6 --header outputs..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_separate_header/callset.json" "GenomicsDB callset.json (separate header index)"
check_file_exists "$meta_temp_dir/genomicsdb_workspace_separate_header/vcfheader.vcf" "GenomicsDB vcfheader.vcf (separate header index)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_separate_header/callset.json" "sample2" "GenomicsDB callset.json sample2 (separate header index)"
check_file_not_contains "$meta_temp_dir/test6_header.log" "no index was provided" "run log (no indexing when an index is given)"
check_file_not_exists "$test_data_dir/separate/header.vcf.gz.tbi" "index next to the header input"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Unindexed .vcf.gz inputs ---
log "Starting TEST 7: .vcf.gz inputs without index files"

mkdir -p "$test_data_dir/no_index"
bgzip -c "$test_data_dir/sample1.g.vcf" > "$test_data_dir/no_index/sample1.g.vcf.gz"
bgzip -c "$test_data_dir/sample2.g.vcf" > "$test_data_dir/no_index/sample2.g.vcf.gz"
bgzip -c "$test_data_dir/sample1.g.vcf" > "$test_data_dir/no_index/header.vcf.gz"

log "Executing $meta_name with unindexed .vcf.gz --variant files..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/sample1.g.vcf.gz" \
  --variant "$test_data_dir/no_index/sample2.g.vcf.gz" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_no_index" 2> "$meta_temp_dir/test7.log"

log "Validating TEST 7 --variant outputs..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_no_index/callset.json" "GenomicsDB callset.json (no indexes)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_no_index/callset.json" "sample1" "GenomicsDB callset.json sample1 (no indexes)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_no_index/callset.json" "sample2" "GenomicsDB callset.json sample2 (no indexes)"
check_file_contains "$meta_temp_dir/test7.log" "no index was provided" "run log (indexing warning)"
check_file_not_exists "$test_data_dir/no_index/sample1.g.vcf.gz.tbi" "index next to the sample1 input"

log "Executing $meta_name with --sample_name_map and an unindexed .vcf.gz --header..."
"$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_no_index_header" \
  --header "$test_data_dir/no_index/header.vcf.gz" 2> "$meta_temp_dir/test7_header.log"

log "Validating TEST 7 --header outputs..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_no_index_header/vcfheader.vcf" "GenomicsDB vcfheader.vcf (no header index)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_no_index_header/callset.json" "sample2" "GenomicsDB callset.json sample2 (no header index)"
check_file_contains "$meta_temp_dir/test7_header.log" "no index was provided" "run log (indexing warning)"
check_file_not_exists "$test_data_dir/no_index/header.vcf.gz.tbi" "index next to the header input"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Index argument errors ---
log "Starting TEST 8: Index argument errors"

log "Executing $meta_name with 2 --variant files and 1 --variant_index file..."
if "$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_mismatch" > "$meta_temp_dir/test8_mismatch.log" 2>&1; then
  log_error "✗ $meta_name did not fail with 2 --variant files and 1 --variant_index file"
  exit 1
fi
check_file_contains "$meta_temp_dir/test8_mismatch.log" "Error: --variant_index must be given once for each file in --variant" "error message"

log "Executing $meta_name with --variant_index and no --variant..."
if "$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_variant_index" > "$meta_temp_dir/test8_variant.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --variant_index and no --variant"
  exit 1
fi
check_file_contains "$meta_temp_dir/test8_variant.log" "Error: --variant_index requires --variant" "error message"

log "Executing $meta_name with --header_index and no --header..."
if "$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --header_index "$test_data_dir/separate_index/header.vcf.gz.tbi" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_header_index" > "$meta_temp_dir/test8_header.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --header_index and no --header"
  exit 1
fi
check_file_contains "$meta_temp_dir/test8_header.log" "Error: --header_index requires --header" "error message"

log "Executing $meta_name with --header and --variant..."
if "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --header "$test_data_dir/sample1.g.vcf" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_header_variant" > "$meta_temp_dir/test8_header_variant.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --header and --variant"
  exit 1
fi
check_file_contains "$meta_temp_dir/test8_header_variant.log" "Error: --header with --variant requires --sample_name" "error message"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Rename staged --variant files with --sample_name ---
log "Starting TEST 9: --sample_name with --variant"

log "Executing $meta_name with --variant, --sample_name and --header..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --sample_name renamed1 \
  --sample_name renamed2 \
  --header "$test_data_dir/sample1.g.vcf" \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_sample_name"

log "Validating TEST 9 outputs..."
check_file_exists "$meta_temp_dir/genomicsdb_workspace_sample_name/callset.json" "GenomicsDB callset.json (sample_name)"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_sample_name/callset.json" "renamed1" "GenomicsDB callset.json renamed1"
check_file_contains "$meta_temp_dir/genomicsdb_workspace_sample_name/callset.json" "renamed2" "GenomicsDB callset.json renamed2"
check_file_not_contains "$meta_temp_dir/genomicsdb_workspace_sample_name/callset.json" "\"sample1\"" "GenomicsDB callset.json (original sample1 name)"

log "Executing $meta_name with fewer --sample_name values than --variant files..."
if "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --sample_name renamed1 \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_sample_name_count" > "$meta_temp_dir/test9_count.log" 2>&1; then
  log_error "✗ $meta_name did not fail with a --sample_name count mismatch"
  exit 1
fi
check_file_contains "$meta_temp_dir/test9_count.log" "Error: --sample_name must be given once for each file in --variant" "error message"

log "Executing $meta_name with --sample_name and --sample_name_map..."
if "$meta_executable" \
  --sample_name_map "$test_data_dir/sample_map.tsv" \
  --sample_name renamed1 \
  --intervals "$test_data_dir/intervals.list" \
  --genomicsdb_workspace_path "$meta_temp_dir/genomicsdb_workspace_sample_name_map" > "$meta_temp_dir/test9_map.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --sample_name and --sample_name_map"
  exit 1
fi
check_file_contains "$meta_temp_dir/test9_map.log" "Error: --sample_name cannot be used together with --sample_name_map" "error message"

log "✅ TEST 9 completed successfully"

print_test_summary "All tests completed successfully"
