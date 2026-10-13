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

mkdir -p "$meta_temp_dir/indexes"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output.vcf" \
  --output_index "$meta_temp_dir/indexes/output.idx"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output.vcf" "output VCF file"
check_file_not_empty "$meta_temp_dir/output.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.vcf" "sample1" "output VCF sample column"
check_file_exists "$meta_temp_dir/indexes/output.idx" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/output.idx" "output index file"
check_file_not_exists "$meta_temp_dir/output.vcf.idx" "index next to the output VCF (moved to --output_index)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: GVCF mode ---
log "Starting TEST 2: Bgzipped GVCF output mode"

log "Executing $meta_name with --emit_ref_confidence GVCF and a .g.vcf.gz output..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output.g.vcf.gz" \
  --output_index "$meta_temp_dir/indexes/output.g.tbi" \
  --emit_ref_confidence GVCF \
  --gvcf_gq_bands 20 \
  --gvcf_gq_bands 60 \
  --floor_blocks

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output.g.vcf.gz" "output GVCF file"
zcat "$meta_temp_dir/output.g.vcf.gz" > "$meta_temp_dir/output.g.vcf"
check_file_not_empty "$meta_temp_dir/output.g.vcf" "output GVCF file"
check_file_contains "$meta_temp_dir/output.g.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.g.vcf" "^##GVCFBlock" "output GVCF block header"
check_file_exists "$meta_temp_dir/indexes/output.g.tbi" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/output.g.tbi" "output index file"
check_file_not_exists "$meta_temp_dir/output.g.vcf.gz.tbi" "index next to the output GVCF (moved to --output_index)"

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
  --output_index "$meta_temp_dir/output_tuned.vcf.idx" \
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

# --- Test Case 5: Bgzipped --dbsnp and --alleles files with separate indexes ---
log "Starting TEST 5: Bgzipped --dbsnp and --alleles files with separate indexes"

create_test_known_sites_vcf "$test_data_dir/sites.vcf" seq1 2000 100
create_test_vcf_gz_with_separate_index "$test_data_dir/sites.vcf" \
  "$test_data_dir/dbsnp.vcf.gz" "$test_data_dir/indexes/dbsnp.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sites.vcf" \
  "$test_data_dir/alleles.vcf.gz" "$test_data_dir/indexes/alleles.vcf.gz.tbi"

log "Executing $meta_name with --dbsnp_index and --alleles_index..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --dbsnp "$test_data_dir/dbsnp.vcf.gz" \
  --dbsnp_index "$test_data_dir/indexes/dbsnp.vcf.gz.tbi" \
  --alleles "$test_data_dir/alleles.vcf.gz" \
  --alleles_index "$test_data_dir/indexes/alleles.vcf.gz.tbi" \
  --output "$meta_temp_dir/output_vcf_inputs.vcf" \
  --output_index "$meta_temp_dir/output_vcf_inputs.vcf.idx"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output_vcf_inputs.vcf" "output VCF file (VCF inputs)"
check_file_contains "$meta_temp_dir/output_vcf_inputs.vcf" "^##fileformat=VCF" "output VCF file header (VCF inputs)"
check_file_contains "$meta_temp_dir/output_vcf_inputs.vcf" "##INFO=<ID=DB," "output VCF DB INFO definition (from --dbsnp)"
check_file_contains "$meta_temp_dir/output_vcf_inputs.vcf" $'^seq1\t100\t' "output VCF record at the force-called site (from --alleles)"
check_file_not_exists "$test_data_dir/dbsnp.vcf.gz.tbi" "index next to the --dbsnp input"
check_file_not_exists "$test_data_dir/alleles.vcf.gz.tbi" "index next to the --alleles input"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Bgzipped --alleles file without an index ---
log "Starting TEST 6: Bgzipped --alleles file without an index"

log "Executing $meta_name with a .vcf.gz --alleles file and no --alleles_index..."
"$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --alleles "$test_data_dir/alleles.vcf.gz" \
  --output "$meta_temp_dir/output_no_index.vcf" \
  --output_index "$meta_temp_dir/output_no_index.vcf.idx" 2> "$meta_temp_dir/test6.log"

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/output_no_index.vcf" "output VCF file (no --alleles index)"
check_file_contains "$meta_temp_dir/output_no_index.vcf" $'^seq1\t100\t' "output VCF record at the force-called site (no --alleles index)"
check_file_contains "$meta_temp_dir/test6.log" "no index was provided" "warning for the missing index"
check_file_not_exists "$test_data_dir/alleles.vcf.gz.tbi" "index next to the --alleles input"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: --dbsnp_index without --dbsnp ---
log "Starting TEST 7: --dbsnp_index without --dbsnp"

log "Executing $meta_name with --dbsnp_index and no --dbsnp (should fail)..."
if "$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --dbsnp_index "$test_data_dir/indexes/dbsnp.vcf.gz.tbi" \
  --output "$meta_temp_dir/output_bad.vcf" > "$meta_temp_dir/test7.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --dbsnp_index and no --dbsnp"
  exit 1
fi
check_file_contains "$meta_temp_dir/test7.log" "Error: --dbsnp_index requires --dbsnp" "error message"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Missing --output_index ---
log "Starting TEST 8: A missing --output_index fails"

if "$meta_executable" \
  --input "$test_data_dir/reads.sorted.bam" \
  --bai "$test_data_dir/reads.sorted.bai" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/no_output_index.vcf" > "$meta_temp_dir/test8.log" 2>&1; then
  log_error "✗ $meta_name did not fail without --output_index"
  exit 1
fi
check_file_contains "$meta_temp_dir/test8.log" "Error: --output_index is required unless --create_output_variant_index is false." "error message"

log "✅ TEST 8 completed successfully"

print_test_summary "All tests"
