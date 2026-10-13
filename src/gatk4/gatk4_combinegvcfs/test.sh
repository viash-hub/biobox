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

# --- Build shared reference ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

# --- Build two per-sample BAMs with distinct read groups/samples ---
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

# --- Generate per-sample GVCFs via a direct HaplotypeCaller invocation ---
log "Running HaplotypeCaller on sample1 to produce a GVCF..."
gatk HaplotypeCaller \
  --input "$test_data_dir/sample1.sorted.bam" \
  --output "$test_data_dir/sample1.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample1.g.vcf" "sample1 GVCF"
check_file_not_empty "$test_data_dir/sample1.g.vcf" "sample1 GVCF"

log "Running HaplotypeCaller on sample2 to produce a GVCF..."
gatk HaplotypeCaller \
  --input "$test_data_dir/sample2.sorted.bam" \
  --output "$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample2.g.vcf" "sample2 GVCF"
check_file_not_empty "$test_data_dir/sample2.g.vcf" "sample2 GVCF"

# --- Test Case 1: Basic combining of two GVCFs ---
log "Starting TEST 1: Combine two per-sample GVCFs"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined.g.vcf" \
  --output_index "$meta_temp_dir/combined.idx"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/combined.g.vcf" "combined GVCF file"
check_file_not_empty "$meta_temp_dir/combined.g.vcf" "combined GVCF file"
check_file_contains "$meta_temp_dir/combined.g.vcf" "^##fileformat=VCF" "combined GVCF header"
check_file_contains "$meta_temp_dir/combined.g.vcf" "sample1" "combined GVCF sample1 column"
check_file_contains "$meta_temp_dir/combined.g.vcf" "sample2" "combined GVCF sample2 column"
check_file_exists "$meta_temp_dir/combined.idx" "output index file (Tribble index for a .vcf)"
check_index_magic "$meta_temp_dir/combined.idx" "TIDX" "output index file"
check_file_not_exists "$meta_temp_dir/combined.g.vcf.idx" "index next to the output GVCF (moved to --output_index)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Additional options ---
log "Starting TEST 2: Combine with additional options"

log "Executing $meta_name with --break_bands_at_multiples_of and --annotation_group..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf;$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_options.g.vcf" \
  --output_index "$meta_temp_dir/combined_options.g.vcf.idx" \
  --break_bands_at_multiples_of 10 \
  --annotation_group StandardAnnotation

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/combined_options.g.vcf" "combined GVCF file with options"
check_file_not_empty "$meta_temp_dir/combined_options.g.vcf" "combined GVCF file with options"
check_file_contains "$meta_temp_dir/combined_options.g.vcf" "^##fileformat=VCF" "combined GVCF (options) header"
check_file_contains "$meta_temp_dir/combined_options.g.vcf" "sample1" "combined GVCF (options) sample1 column"
check_file_contains "$meta_temp_dir/combined_options.g.vcf" "sample2" "combined GVCF (options) sample2 column"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Annotation and output-shaping options ---
log "Starting TEST 3: Combine with annotation and output-shaping options"

log "Executing $meta_name with --annotation, --annotations_to_exclude, --sites_only_vcf_output, --variant_output_filtering, and --ignore_variants_starting_outside_interval..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_annotations.g.vcf" \
  --output_index "$meta_temp_dir/combined_annotations.g.vcf.idx" \
  --annotation ChromosomeCounts \
  --annotations_to_exclude InbreedingCoeff \
  --sites_only_vcf_output \
  --variant_output_filtering ANYWHERE \
  --ignore_variants_starting_outside_interval

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/combined_annotations.g.vcf" "combined GVCF file with annotation options"
check_file_not_empty "$meta_temp_dir/combined_annotations.g.vcf" "combined GVCF file with annotation options"
check_file_contains "$meta_temp_dir/combined_annotations.g.vcf" "^##fileformat=VCF" "combined GVCF (annotations) header"
check_file_contains "$meta_temp_dir/combined_annotations.g.vcf" "##INFO=<ID=AC," "combined GVCF (annotations) AC INFO definition"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/combined_annotations.g.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* || "$chrom_line" == *"sample2"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Shared GATK engine options ---
log "Starting TEST 4: Shared GATK engine options"

log "Executing $meta_name with --interval_padding and --create_output_variant_index false..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_engine.g.vcf" \
  --interval_padding 10 \
  --create_output_variant_index false

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/combined_engine.g.vcf" "combined GVCF file (engine options)"
check_file_not_empty "$meta_temp_dir/combined_engine.g.vcf" "combined GVCF file (engine options)"
check_file_not_exists "$meta_temp_dir/combined_engine.g.vcf.idx" "combined GVCF index (should not exist with --create_output_variant_index false)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Index files in a different directory ---
log "Starting TEST 5: --variant_index and --dbsnp_index in a different directory"

mkdir -p "$test_data_dir/separate"

create_test_vcf_gz_with_separate_index "$test_data_dir/sample1.g.vcf" \
  "$test_data_dir/separate/sample1.g.vcf.gz" "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sample2.g.vcf" \
  "$test_data_dir/separate/sample2.g.vcf.gz" "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi"
create_test_known_sites_vcf "$test_data_dir/dbsnp.vcf" seq1 2000 100
create_test_vcf_gz_with_separate_index "$test_data_dir/dbsnp.vcf" \
  "$test_data_dir/separate/dbsnp.vcf.gz" "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi"

log "Executing $meta_name with .vcf.gz inputs and separate .tbi index files..."
"$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --variant_index "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi" \
  --dbsnp "$test_data_dir/separate/dbsnp.vcf.gz" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_separate.g.vcf" \
  --output_index "$meta_temp_dir/combined_separate.g.vcf.idx" 2> "$meta_temp_dir/test5.log"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/combined_separate.g.vcf" "combined GVCF file (separate indexes)"
check_file_contains "$meta_temp_dir/combined_separate.g.vcf" "^##fileformat=VCF" "combined GVCF (separate indexes) header"
check_file_contains "$meta_temp_dir/combined_separate.g.vcf" "sample1" "combined GVCF (separate indexes) sample1 column"
check_file_contains "$meta_temp_dir/combined_separate.g.vcf" "sample2" "combined GVCF (separate indexes) sample2 column"
check_file_not_contains "$meta_temp_dir/test5.log" "no index was provided" "run log (no indexing when indexes are given)"
check_file_not_exists "$test_data_dir/separate/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/separate/sample2.g.vcf.gz.tbi" "index next to the sample2 input"
check_file_not_exists "$test_data_dir/separate/dbsnp.vcf.gz.tbi" "index next to the dbsnp input"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Unindexed .vcf.gz inputs ---
log "Starting TEST 6: .vcf.gz inputs without index files"

mkdir -p "$test_data_dir/no_index"
bgzip -c "$test_data_dir/sample1.g.vcf" > "$test_data_dir/no_index/sample1.g.vcf.gz"
bgzip -c "$test_data_dir/sample2.g.vcf" > "$test_data_dir/no_index/sample2.g.vcf.gz"
bgzip -c "$test_data_dir/dbsnp.vcf" > "$test_data_dir/no_index/dbsnp.vcf.gz"

log "Executing $meta_name with unindexed .vcf.gz inputs..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/sample1.g.vcf.gz" \
  --variant "$test_data_dir/no_index/sample2.g.vcf.gz" \
  --dbsnp "$test_data_dir/no_index/dbsnp.vcf.gz" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_no_index.g.vcf" \
  --output_index "$meta_temp_dir/combined_no_index.g.vcf.idx" 2> "$meta_temp_dir/test6.log"

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/combined_no_index.g.vcf" "combined GVCF file (no indexes)"
check_file_contains "$meta_temp_dir/combined_no_index.g.vcf" "sample1" "combined GVCF (no indexes) sample1 column"
check_file_contains "$meta_temp_dir/combined_no_index.g.vcf" "sample2" "combined GVCF (no indexes) sample2 column"
check_file_contains "$meta_temp_dir/test6.log" "no index was provided" "run log (indexing warning)"
check_file_not_exists "$test_data_dir/no_index/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/no_index/dbsnp.vcf.gz.tbi" "index next to the dbsnp input"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Index argument errors ---
log "Starting TEST 7: Index argument errors"

log "Executing $meta_name with 2 --variant files and 1 --variant_index file..."
if "$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_mismatch.g.vcf" \
  --output_index "$meta_temp_dir/combined_mismatch.g.vcf.idx" > "$meta_temp_dir/test7_mismatch.log" 2>&1; then
  log_error "✗ $meta_name did not fail with 2 --variant files and 1 --variant_index file"
  exit 1
fi
check_file_contains "$meta_temp_dir/test7_mismatch.log" "Error: --variant_index must be given once for each file in --variant" "error message"

log "Executing $meta_name with --dbsnp_index and no --dbsnp..."
if "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/combined_dbsnp_index.g.vcf" \
  --output_index "$meta_temp_dir/combined_dbsnp_index.g.vcf.idx" > "$meta_temp_dir/test7_dbsnp.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --dbsnp_index and no --dbsnp"
  exit 1
fi
check_file_contains "$meta_temp_dir/test7_dbsnp.log" "Error: --dbsnp_index requires --dbsnp" "error message"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Bgzipped output with an index in a different directory ---
log "Starting TEST 8: .g.vcf.gz output with --output_index in a different directory"

mkdir -p "$meta_temp_dir/gz_output" "$meta_temp_dir/gz_output_index"

log "Executing $meta_name with a .g.vcf.gz --output and --output_index in a different directory..."
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/gz_output/combined.g.vcf.gz" \
  --output_index "$meta_temp_dir/gz_output_index/combined.g.vcf.gz.tbi"

log "Validating TEST 8 outputs..."
check_file_exists "$meta_temp_dir/gz_output/combined.g.vcf.gz" "combined GVCF file (bgzipped)"
zcat "$meta_temp_dir/gz_output/combined.g.vcf.gz" > "$meta_temp_dir/gz_output_index/combined.g.vcf"
check_file_contains "$meta_temp_dir/gz_output_index/combined.g.vcf" "sample1" "combined GVCF (bgzipped) sample1 column"
check_file_contains "$meta_temp_dir/gz_output_index/combined.g.vcf" "sample2" "combined GVCF (bgzipped) sample2 column"
check_file_exists "$meta_temp_dir/gz_output_index/combined.g.vcf.gz.tbi" "output index file (tabix index for a .vcf.gz)"
check_index_magic "$meta_temp_dir/gz_output_index/combined.g.vcf.gz.tbi" "TBI" "output index file"

log "Checking that no index was left next to the output GVCF..."
found_files=$(find "$meta_temp_dir/gz_output" -type f -exec basename {} \; | sort)
if [[ "$found_files" == "combined.g.vcf.gz" ]]; then
  log "✓ Output directory contains only the output GVCF"
else
  log_error "✗ Output directory contains unexpected files: $(echo "$found_files" | tr '\n' ' ')"
  exit 1
fi

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Missing --output_index ---
log "Starting TEST 9: A missing --output_index fails"

if "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/no_output_index.g.vcf" > "$meta_temp_dir/test9.log" 2>&1; then
  log_error "✗ $meta_name did not fail without --output_index"
  exit 1
fi
check_file_contains "$meta_temp_dir/test9.log" "Error: --output_index is required unless --create_output_variant_index is false." "error message"
check_file_not_exists "$meta_temp_dir/no_output_index.g.vcf" "output GVCF (no run without --output_index)"

log "✅ TEST 9 completed successfully"

print_test_summary "All tests"
