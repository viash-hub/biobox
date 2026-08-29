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

# --- Build shared reference fixtures ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

# --- Build a GVCF fixture by running HaplotypeCaller directly ---
create_test_sam_reads "$test_data_dir/reads.sam" seq1 2000 rg1 sample1
check_file_exists "$test_data_dir/reads.sam" "hand-crafted SAM reads"

sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bam" "sorted BAM file"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"

log "Running HaplotypeCaller in GVCF mode to build the input fixture..."
gatk HaplotypeCaller \
  --input "$test_data_dir/reads.sorted.bam" \
  --output "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample.g.vcf" "input GVCF fixture"
check_file_contains "$test_data_dir/sample.g.vcf" "<NON_REF>" "input GVCF fixture NON_REF placeholder allele"

# --- Test Case 1: Default parameters ---
log "Starting TEST 1: Default parameters"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output.vcf"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output.vcf" "output VCF file"
check_file_not_empty "$meta_temp_dir/output.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.vcf" "sample1" "output VCF sample column"
check_file_not_contains "$meta_temp_dir/output.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Options ---
log "Starting TEST 2: Additional options"

log "Executing $meta_name with additional options..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output_options.vcf" \
  --include_non_variant_sites \
  --stand_call_conf 10.0 \
  --sample_ploidy 2

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output_options.vcf" "output VCF file with options"
check_file_not_empty "$meta_temp_dir/output_options.vcf" "output VCF file with options"
check_file_contains "$meta_temp_dir/output_options.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_not_contains "$meta_temp_dir/output_options.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Annotation and output-shaping options ---
log "Starting TEST 3: Additional annotation and output-shaping options"

log "Executing $meta_name with --annotation, --annotations_to_exclude, --call_genotypes, --sites_only_vcf_output, and --variant_output_filtering..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output_annotations.vcf" \
  --annotation ChromosomeCounts \
  --annotations_to_exclude InbreedingCoeff \
  --call_genotypes \
  --sites_only_vcf_output \
  --variant_output_filtering ANYWHERE

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output_annotations.vcf" "output VCF file with annotation options"
check_file_not_empty "$meta_temp_dir/output_annotations.vcf" "output VCF file with annotation options"
check_file_contains "$meta_temp_dir/output_annotations.vcf" "^##fileformat=VCF" "output VCF (annotations) header"
check_file_contains "$meta_temp_dir/output_annotations.vcf" "##INFO=<ID=AC," "output VCF (annotations) AC INFO definition"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/output_annotations.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
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
  --variant "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output_engine.vcf" \
  --interval_padding 10 \
  --create_output_variant_index false

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_empty "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_exists "$meta_temp_dir/output_engine.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Test GenomicsDB workspace input ---
log "Starting TEST 5: GenomicsDB workspace input"

log "Creating a GenomicsDB workspace from the input GVCF..."
printf 'sample1\t%s\n' "$test_data_dir/sample.g.vcf" > "$test_data_dir/sample_map.tsv"
gatk GenomicsDBImport \
  --genomicsdb-workspace-path "$meta_temp_dir/genomicsdb_workspace" \
  --batch-size 50 \
  --sample-name-map "$test_data_dir/sample_map.tsv" \
  --reader-threads 1 \
  --tmp-dir "$meta_temp_dir" \
  --intervals "seq1:1-2000" \
  --reference "$test_data_dir/reference.fasta"

log "Executing $meta_name with GenomicsDB workspace input..."
"$meta_executable" \
  --variant "$meta_temp_dir/genomicsdb_workspace" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --output "$meta_temp_dir/output_genomicsdb.vcf"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output_genomicsdb.vcf" "output VCF file (GenomicsDB input)"
check_file_not_empty "$meta_temp_dir/output_genomicsdb.vcf" "output VCF file (GenomicsDB input)"
check_file_contains "$meta_temp_dir/output_genomicsdb.vcf" "^##fileformat=VCF" "output VCF file header (GenomicsDB input)"
check_file_contains "$meta_temp_dir/output_genomicsdb.vcf" "sample1" "output VCF sample column (GenomicsDB input)"
check_file_not_contains "$meta_temp_dir/output_genomicsdb.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele (GenomicsDB input)"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests completed successfully"
