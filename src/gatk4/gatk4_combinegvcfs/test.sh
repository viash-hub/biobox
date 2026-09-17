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
  --output "$meta_temp_dir/combined.g.vcf"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/combined.g.vcf" "combined GVCF file"
check_file_not_empty "$meta_temp_dir/combined.g.vcf" "combined GVCF file"
check_file_contains "$meta_temp_dir/combined.g.vcf" "^##fileformat=VCF" "combined GVCF header"
check_file_contains "$meta_temp_dir/combined.g.vcf" "sample1" "combined GVCF sample1 column"
check_file_contains "$meta_temp_dir/combined.g.vcf" "sample2" "combined GVCF sample2 column"

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

print_test_summary "All tests"
