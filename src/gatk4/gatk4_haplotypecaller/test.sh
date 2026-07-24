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

# --- Build shared reference and BAM fixtures ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

create_test_sam_reads "$test_data_dir/reads.sam" seq1 2000 rg1 sample1
check_file_exists "$test_data_dir/reads.sam" "hand-crafted SAM reads"

sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bam" "sorted BAM file"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"

# --- Test Case 1: Default parameters ---
log "Starting TEST 1: Default parameters"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output.vcf"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output.vcf" "output VCF file"
check_file_not_empty "$meta_temp_dir/output.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.vcf" "sample1" "output VCF sample column"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: GVCF mode ---
log "Starting TEST 2: GVCF output mode"

log "Executing $meta_name with --emit_ref_confidence GVCF..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output.g.vcf" \
  --emit_ref_confidence GVCF \
  --gvcf_gq_bands 20 \
  --gvcf_gq_bands 60 \
  --floor_blocks

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output.g.vcf" "output GVCF file"
check_file_not_empty "$meta_temp_dir/output.g.vcf" "output GVCF file"
check_file_contains "$meta_temp_dir/output.g.vcf" "^##fileformat=VCF" "output GVCF file header"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Assembly/genotyping tuning options ---
log "Starting TEST 3: Assembly and genotyping tuning options"

log "Executing $meta_name with sample-ploidy, output-mode, annotation and other tuning options..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output_tuned.vcf" \
  --sample_ploidy 2 \
  --output_mode EMIT_ALL_CONFIDENT_SITES \
  --dont_use_soft_clipped_bases \
  --minimum_mapping_quality 20 \
  --base_quality_score_threshold 18 \
  --assembly_region_padding 100 \
  --kmer_size 10 \
  --kmer_size 25 \
  --min_pruning 2 \
  --active_probability_threshold 0.002 \
  --annotation ChromosomeCounts \
  --annotations_to_exclude InbreedingCoeff

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output_tuned.vcf" "output VCF file"
check_file_not_empty "$meta_temp_dir/output_tuned.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output_tuned.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output_tuned.vcf" "sample1" "output VCF sample column"
check_file_contains "$meta_temp_dir/output_tuned.vcf" "##INFO=<ID=AC," "output VCF AC INFO definition (from --annotation ChromosomeCounts)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Shared GATK engine options ---
log "Starting TEST 4: Shared GATK engine options"

log "Executing $meta_name with --interval_padding, --sites_only_vcf_output, --create_output_variant_index false..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output_engine.vcf" \
  --interval_padding 10 \
  --sites_only_vcf_output \
  --create_output_variant_index false

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_empty "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_exists "$meta_temp_dir/output_engine.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/output_engine.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
