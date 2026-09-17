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

print_test_summary "All tests completed successfully"
