#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---

# ATCG-repeat reference (seq1, 2000 bp) with 6 reads, for the tests with
# known-site files and intervals
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"
reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

create_test_sam_reads "$test_data_dir/reads.sam" seq1 2000 rg1 sample1
sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"
atcg_args=(
  --input "$test_data_dir/reads.sorted.bam"
  --read_index "$test_data_dir/reads.sorted.bai"
  "${reference_args[@]}"
)

# Somatic reference (chr1, 500 bp) with a SNV (chr1:250 T>A) in half of the
# reads of snv_sample, all reads of normal_sample are reference. The
# snv_sample reads are also written as a CRAM with the .crai in another
# directory.
mkdir -p "$test_data_dir/snv"
create_test_somatic_reference "$test_data_dir/snv/reference.fasta"
create_test_somatic_bam "$test_data_dir/snv/sample.bam" snv_sample rg1 mixed
create_test_somatic_bam "$test_data_dir/snv/normal.bam" normal_sample rg2 ref
create_test_cram "$test_data_dir/snv/sample.bam" "$test_data_dir/snv/reference.fasta" \
  "$test_data_dir/cram/sample.cram" "$test_data_dir/cram_index/sample_index.crai"
snv_reference_args=(
  --reference "$test_data_dir/snv/reference.fasta"
  --reference_fai "$test_data_dir/snv/reference.fasta.fai"
  --reference_dict "$test_data_dir/snv/reference.dict"
)
snv_args=(
  --input "$test_data_dir/snv/sample.bam"
  --read_index "$test_data_dir/snv/sample.bai"
  "${snv_reference_args[@]}"
)

# Bgzipped copies of a one-site VCF (seq1:100) for --dbsnp and --alleles,
# with their indexes in another directory
create_test_known_sites_vcf --no-index "$test_data_dir/sites.vcf" seq1 2000 100
create_test_vcf_gz_with_separate_index "$test_data_dir/sites.vcf" \
  "$test_data_dir/dbsnp.vcf.gz" "$test_data_dir/indexes/dbsnp.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/sites.vcf" \
  "$test_data_dir/alleles.vcf.gz" "$test_data_dir/indexes/alleles.vcf.gz.tbi"

# Bgzipped population callset with the SNV, with its index in another
# directory
create_test_known_sites_vcf --no-index "$test_data_dir/snv/population.vcf" chr1 500 250 T A
create_test_vcf_gz_with_separate_index "$test_data_dir/snv/population.vcf" \
  "$test_data_dir/snv/population.vcf.gz" "$test_data_dir/snv/population_index/population.vcf.gz.tbi"

# --- Test Case 1: Basic run (defaults) ---
log "Starting TEST 1: Basic run (defaults)"

mkdir -p "$meta_temp_dir/indexes"

"$meta_executable" \
  "${snv_args[@]}" \
  --output "$meta_temp_dir/output.vcf" \
  --output_index "$meta_temp_dir/indexes/output.idx"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/output.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.vcf" $'^#CHROM.*\tsnv_sample$' "output VCF sample column"
check_vcf_has_record "$meta_temp_dir/output.vcf" $'^chr1\t250\t.\tT\tA\t' "SNV chr1:250 T>A"
check_index_moved "$meta_temp_dir/indexes/output.idx" "$meta_temp_dir/output.vcf.idx" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/output.idx" "output index file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Inputs and indexes"

# Run A: CRAM input with its index in another directory
log "Executing $meta_name with a CRAM --input and a --read_index in another directory..."
"$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample_index.crai" \
  "${snv_reference_args[@]}" \
  --output "$meta_temp_dir/snv_cram.vcf" \
  --output_index "$meta_temp_dir/snv_cram.vcf.idx"

check_vcf_has_record "$meta_temp_dir/snv_cram.vcf" $'^chr1\t250\t.\tT\tA\t' "SNV chr1:250 T>A called from the CRAM input"
bam_records=$(grep -v '^#' "$meta_temp_dir/output.vcf")
cram_records=$(grep -v '^#' "$meta_temp_dir/snv_cram.vcf")
if [[ "$cram_records" == "$bam_records" ]]; then
  log "✓ CRAM input gives the same calls as the BAM input of TEST 1"
else
  log_error "✗ CRAM input gives different calls from the BAM input of TEST 1"
  exit 1
fi

# Run B: two BAM inputs from different samples. --stand_call_conf is the
# Viash alternative name of --standard_min_confidence_threshold_for_calling.
log "Executing $meta_name with two --input and two --read_index files..."
"$meta_executable" \
  --input "$test_data_dir/snv/sample.bam" \
  --input "$test_data_dir/snv/normal.bam" \
  --read_index "$test_data_dir/snv/sample.bai" \
  --read_index "$test_data_dir/snv/normal.bai" \
  "${snv_reference_args[@]}" \
  --output "$meta_temp_dir/two_samples.vcf" \
  --output_index "$meta_temp_dir/two_samples.vcf.idx" \
  --stand_call_conf 20

check_file_contains "$meta_temp_dir/two_samples.vcf" $'^#CHROM.*\tFORMAT\tnormal_sample\tsnv_sample$' "output VCF columns for both samples"
check_vcf_has_record "$meta_temp_dir/two_samples.vcf" $'^chr1\t250\t.\tT\tA\t.*\t0/0:[^\t]*\t0/1:' "SNV chr1:250 T>A with a hom-ref and a het genotype"
check_file_contains "$meta_temp_dir/two_samples.vcf" " --standard-min-confidence-threshold-for-calling 20.0" "GATK command line (from the --stand_call_conf alternative)"

# Run C: bgzipped --dbsnp and --alleles files with separate indexes
log "Executing $meta_name with --dbsnp_index and --alleles_index..."
"$meta_executable" \
  "${atcg_args[@]}" \
  --dbsnp "$test_data_dir/dbsnp.vcf.gz" \
  --dbsnp_index "$test_data_dir/indexes/dbsnp.vcf.gz.tbi" \
  --alleles "$test_data_dir/alleles.vcf.gz" \
  --alleles_index "$test_data_dir/indexes/alleles.vcf.gz.tbi" \
  --output "$meta_temp_dir/output_vcf_inputs.vcf" \
  --output_index "$meta_temp_dir/output_vcf_inputs.vcf.idx"

check_file_contains "$meta_temp_dir/output_vcf_inputs.vcf" "^##fileformat=VCF" "output VCF file header (VCF inputs)"
check_file_contains "$meta_temp_dir/output_vcf_inputs.vcf" "##INFO=<ID=DB," "output VCF DB INFO definition (from --dbsnp)"
check_vcf_has_record "$meta_temp_dir/output_vcf_inputs.vcf" $'^seq1\t100\t' "output VCF record at the force-called site (from --alleles)"
check_file_not_exists "$test_data_dir/dbsnp.vcf.gz.tbi" "index next to the --dbsnp input"
check_file_not_exists "$test_data_dir/alleles.vcf.gz.tbi" "index next to the --alleles input"

# Run D: bgzipped --alleles file without an index, with the engine options
# and a --bam_output path without --create_bam_output
log "Executing $meta_name with an unindexed .vcf.gz --alleles file, engine options and --bam_output without --create_bam_output..."
"$meta_executable" \
  "${atcg_args[@]}" \
  --alleles "$test_data_dir/alleles.vcf.gz" \
  --output "$meta_temp_dir/output_no_index.vcf" \
  --interval_padding 10 \
  --sites_only_vcf_output \
  --create_output_variant_index false \
  --bam_output "$meta_temp_dir/no_bamout/haplotypes.bam" \
  --bam_output_index "$meta_temp_dir/no_bamout/haplotypes.bai" \
  2> "$meta_temp_dir/no_index.log"

check_file_not_empty "$meta_temp_dir/output_no_index.vcf" "output VCF file (no --alleles index)"
check_vcf_has_record "$meta_temp_dir/output_no_index.vcf" $'^seq1\t100\t' "output VCF record at the force-called site (no --alleles index)"
check_file_contains "$meta_temp_dir/no_index.log" "no index was provided" "warning for the missing index"
check_file_not_exists "$test_data_dir/alleles.vcf.gz.tbi" "index next to the --alleles input"
check_file_not_exists "$meta_temp_dir/output_no_index.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/output_no_index.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (--sites_only_vcf_output): $chrom_line"
fi
check_file_not_exists "$meta_temp_dir/no_bamout/haplotypes.bam" "output BAM file (no --create_bam_output)"
check_file_not_exists "$meta_temp_dir/no_bamout/haplotypes.bai" "output BAM index file (no --create_bam_output)"

# Run E: population callset with a separate index
log "Executing $meta_name with --population_callset and --population_callset_index..."
"$meta_executable" \
  "${snv_args[@]}" \
  --output "$meta_temp_dir/population.vcf" \
  --output_index "$meta_temp_dir/population.vcf.idx" \
  --population_callset "$test_data_dir/snv/population.vcf.gz" \
  --population_callset_index "$test_data_dir/snv/population_index/population.vcf.gz.tbi" \
  --num_reference_samples_if_no_call 2 \
  --use_posteriors_to_calculate_qual

check_file_contains "$meta_temp_dir/population.vcf" "##INFO=<ID=PG," "output VCF genotype prior INFO definition (from --population_callset)"
check_file_contains "$meta_temp_dir/population.vcf" " --use-posteriors-to-calculate-qual true" "GATK command line (from --use_posteriors_to_calculate_qual)"
check_vcf_has_record "$meta_temp_dir/population.vcf" $'^chr1\t250\t.\tT\tA\t' "SNV chr1:250 T>A"
check_file_not_exists "$test_data_dir/snv/population.vcf.gz.tbi" "index next to the --population_callset input"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Multiple --intervals files (BED and bgzipped VCF) with --exclude_intervals"

# The territory is the union of seq1:1-300 (BED) and seq1:1550 (VCF), minus
# seq1:101-150 (BED). In GVCF mode, the records cover all of the territory.
printf 'seq1\t0\t300\n' > "$test_data_dir/intervals_a.bed"
printf 'seq1\t100\t150\n' > "$test_data_dir/exclude.bed"
create_test_known_sites_vcf --no-index "$test_data_dir/intervals_site.vcf" seq1 2000 1550
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

"$meta_executable" \
  "${atcg_args[@]}" \
  --output "$meta_temp_dir/output_intervals.g.vcf" \
  --output_index "$meta_temp_dir/indexes/output_intervals.g.vcf.idx" \
  --emit_ref_confidence GVCF \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule UNION \
  --exclude_intervals "$test_data_dir/exclude.bed"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output_intervals.g.vcf" "output GVCF (multiple intervals)"
check_file_exists "$meta_temp_dir/indexes/output_intervals.g.vcf.idx" "output index (multiple intervals)"

# Merge the abutting GVCF records (POS to INFO/END) into regions
covered=$(grep -v '^#' "$meta_temp_dir/output_intervals.g.vcf" | awk -F'\t' '
  {
    start = $2; end = $2
    if (match($8, /(^|;)END=[0-9]+/)) { end = substr($8, RSTART, RLENGTH); sub(/.*END=/, "", end) }
    if (n > 0 && start == region_end + 1) { region_end = end; next }
    if (n > 0) printf "%d-%d ", region_start, region_end
    region_start = start; region_end = end; n++
  }
  END { if (n > 0) printf "%d-%d ", region_start, region_end }')
expected_covered="1-100 151-300 1550-1550 "
if [[ "$covered" == "$expected_covered" ]]; then
  log "✓ output GVCF only covers the selected territory: $covered"
else
  log_error "✗ output GVCF covers '$covered', expected '$expected_covered'"
  exit 1
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Outputs ---
log "Starting TEST 4: Bgzipped GVCF output and haplotype BAM output"

"$meta_executable" \
  "${snv_args[@]}" \
  --output "$meta_temp_dir/output.g.vcf.gz" \
  --output_index "$meta_temp_dir/indexes/output.g.tbi" \
  --emit_ref_confidence GVCF \
  --gvcf_gq_bands 20 \
  --gvcf_gq_bands 60 \
  --floor_blocks \
  --create_bam_output \
  --bam_output "$meta_temp_dir/bamout/haplotypes.bam" \
  --bam_output_index "$meta_temp_dir/bamout_index/haplotypes_index.bai" \
  --bam_writer_type ALL_POSSIBLE_HAPLOTYPES

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output.g.vcf.gz" "output GVCF file"
zcat "$meta_temp_dir/output.g.vcf.gz" > "$meta_temp_dir/output.g.vcf"
check_file_not_empty "$meta_temp_dir/output.g.vcf" "output GVCF file"
check_file_contains "$meta_temp_dir/output.g.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.g.vcf" "^##GVCFBlock20-60=" "output GVCF block header (from --gvcf_gq_bands)"
check_file_contains "$meta_temp_dir/output.g.vcf" " --floor-blocks true" "GATK command line (from --floor_blocks)"
check_index_moved "$meta_temp_dir/indexes/output.g.tbi" "$meta_temp_dir/output.g.vcf.gz.tbi" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/output.g.tbi" "output index file"

check_file_not_empty "$meta_temp_dir/bamout/haplotypes.bam" "output BAM file"
check_index_moved "$meta_temp_dir/bamout_index/haplotypes_index.bai" "$meta_temp_dir/bamout/haplotypes.bai" "output BAM index file"
check_file_not_empty "$meta_temp_dir/bamout_index/haplotypes_index.bai" "output BAM index file"
if samtools view -H "$meta_temp_dir/bamout/haplotypes.bam" | grep -q $'^@RG\tID:ArtificialHaplotypeRG'; then
  log "✓ output BAM has the artificial haplotype read group"
else
  log_error "✗ output BAM has no artificial haplotype read group"
  exit 1
fi

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Tool options ---
log "Starting TEST 5: Assembly, genotyping, mapping quality, phasing and annotation options"

printf 'snv_sample\t0.0\n' > "$test_data_dir/snv/contamination.tsv"

# Run A: the options do not remove the SNV. --annotation ChromosomeCounts
# adds AC, --disable_tool_default_annotations removes QD.
log "Executing $meta_name with the tool options..."
"$meta_executable" \
  "${snv_args[@]}" \
  --output "$meta_temp_dir/calling_options.vcf" \
  --output_index "$meta_temp_dir/calling_options.vcf.idx" \
  --sample_ploidy 2 \
  --output_mode EMIT_ALL_CONFIDENT_SITES \
  --dont_use_soft_clipped_bases \
  --base_quality_score_threshold 18 \
  --assembly_region_padding 100 \
  --kmer_size 10 \
  --kmer_size 25 \
  --min_pruning 2 \
  --active_probability_threshold 0.002 \
  --annotation ChromosomeCounts \
  --annotation PossibleDeNovo \
  --annotations_to_exclude InbreedingCoeff \
  --minimum_mapping_quality 30 \
  --mapping_quality_threshold_for_genotyping 30 \
  --maximum_mapping_quality 60 \
  --do_not_run_physical_phasing \
  --max_mnp_distance 1 \
  --max_genotype_count 512 \
  --contamination_fraction_per_sample_file "$test_data_dir/snv/contamination.tsv" \
  --denovo_depth_threshold 5 \
  --denovo_parent_gq_threshold 10 \
  --disable_tool_default_annotations \
  --standard_min_confidence_threshold_for_calling 10

check_file_contains "$meta_temp_dir/calling_options.vcf" $'^#CHROM.*\tsnv_sample$' "output VCF sample column"
check_vcf_has_record "$meta_temp_dir/calling_options.vcf" $'^chr1\t250\t.\tT\tA\t' "SNV chr1:250 T>A"
check_file_contains "$meta_temp_dir/calling_options.vcf" "##INFO=<ID=AC," "output VCF AC INFO definition (from --annotation ChromosomeCounts)"
check_file_not_contains "$meta_temp_dir/calling_options.vcf" "##INFO=<ID=QD," "output VCF QD INFO definition (removed by --disable_tool_default_annotations)"
check_file_contains "$meta_temp_dir/calling_options.vcf" " --do-not-run-physical-phasing true" "GATK command line (from --do_not_run_physical_phasing)"
check_file_contains "$meta_temp_dir/calling_options.vcf" " --mapping-quality-threshold-for-genotyping 30" "GATK command line (from --mapping_quality_threshold_for_genotyping)"
check_file_contains "$meta_temp_dir/calling_options.vcf" " --max-mnp-distance 1" "GATK command line (from --max_mnp_distance)"
check_file_contains "$meta_temp_dir/calling_options.vcf" " --standard-min-confidence-threshold-for-calling 10.0" "GATK command line (from --standard_min_confidence_threshold_for_calling)"

# Run B: all reads have mapping quality 60, so --maximum_mapping_quality 59
# removes all of them and the SNV of TEST 1 is not called
log "Executing $meta_name with --maximum_mapping_quality below the read mapping quality..."
"$meta_executable" \
  "${snv_args[@]}" \
  --output "$meta_temp_dir/max_mapq.vcf" \
  --output_index "$meta_temp_dir/max_mapq.vcf.idx" \
  --maximum_mapping_quality 59
check_vcf_lacks_record "$meta_temp_dir/max_mapq.vcf" $'^chr1\t250\t' "SNV chr1:250 (all reads removed by --maximum_mapping_quality)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Argument errors ---
log "Starting TEST 6: Argument errors"

# The command for the runs that must fail
fail_cmd=(
  "$meta_executable"
  "${snv_args[@]}"
  --output "$meta_temp_dir/fail.vcf"
  --output_index "$meta_temp_dir/fail.vcf.idx"
)

log "Executing $meta_name with --dbsnp_index and no --dbsnp (should fail)..."
check_fails_with "Error: --dbsnp_index requires --dbsnp" "$meta_temp_dir/fail_dbsnp_index.log" \
  "${fail_cmd[@]}" --dbsnp_index "$test_data_dir/indexes/dbsnp.vcf.gz.tbi"

log "Executing $meta_name without --output_index (should fail)..."
check_fails_with "Error: --output_index is required unless --create_output_variant_index is false." "$meta_temp_dir/fail_output_index.log" \
  "$meta_executable" \
  "${snv_args[@]}" \
  --output "$meta_temp_dir/no_output_index.vcf"

log "Executing $meta_name with two --input files and one --read_index file (should fail)..."
check_fails_with "Error: --read_index must be given once for each file in --input" "$meta_temp_dir/fail_read_index.log" \
  "$meta_executable" \
  --input "$test_data_dir/snv/sample.bam" \
  --input "$test_data_dir/snv/normal.bam" \
  --read_index "$test_data_dir/snv/sample.bai" \
  "${snv_reference_args[@]}" \
  --output "$meta_temp_dir/fail.vcf" \
  --output_index "$meta_temp_dir/fail.vcf.idx"

log "Executing $meta_name with invalid --bam_output arguments (should fail)..."
check_fails_with "Error: --create_bam_output requires --bam_output." "$meta_temp_dir/fail_bam_output.log" \
  "${fail_cmd[@]}" --create_bam_output
check_fails_with "Error: --bam_output_index is required unless --create_output_bam_index is false." "$meta_temp_dir/fail_bam_output_index.log" \
  "${fail_cmd[@]}" --create_bam_output --bam_output "$meta_temp_dir/fail.bam"
check_fails_with "Error: --bam_output must be a BAM file, GATK cannot write CRAM here." "$meta_temp_dir/fail_bam_output_cram.log" \
  "${fail_cmd[@]}" --create_bam_output --bam_output "$meta_temp_dir/fail.cram" --bam_output_index "$meta_temp_dir/fail.crai"
check_file_not_contains "$meta_temp_dir/fail_bam_output_cram.log" "HaplotypeCaller done" "run log (GATK did not run)"

log "Executing $meta_name with --population_callset_index and no --population_callset (should fail)..."
check_fails_with "Error: --population_callset_index requires --population_callset." "$meta_temp_dir/fail_population_index.log" \
  "${fail_cmd[@]}" --population_callset_index "$test_data_dir/snv/population_index/population.vcf.gz.tbi"

check_file_not_exists "$meta_temp_dir/fail.vcf" "output VCF of the runs that must fail"
check_file_not_exists "$meta_temp_dir/no_output_index.vcf" "output VCF without --output_index"

log "✅ TEST 6 completed successfully"

print_test_summary "All tests"
