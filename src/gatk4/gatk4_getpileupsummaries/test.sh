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

# --- Build shared reference, VCF and BAM fixtures ---
# A 200bp single-contig reference made entirely of "A" bases, so any read
# with a different base at a site unambiguously carries the alternate allele.
# Two sites are "common SNPs": position 50 (A>G, population AF 0.3) and
# position 150 (A>T, population AF 0.4).
log "Writing all-A reference FASTA..."
a80=$(printf 'A%.0s' {1..80})
a40=$(printf 'A%.0s' {1..40})
printf '>chr1\n%s\n%s\n%s\n' "$a80" "$a80" "$a40" > "$test_data_dir/reference.fasta"
gatk CreateSequenceDictionary -R "$test_data_dir/reference.fasta" -O "$test_data_dir/reference.dict" --VERBOSITY ERROR
create_test_fasta_fai "$test_data_dir/reference.fasta" "$test_data_dir/reference.fasta.fai"

log "Writing common SNPs VCF..."
cat > "$test_data_dir/common_snps.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=200>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	50	.	A	G	.	.	AF=0.3
chr1	150	.	A	T	.	.	AF=0.4
VCF

# Reads 1-3 start at position 1 and cover position 50 (offset 49): 2 carry
# the reference "A" and 1 carries the alternate "G". Reads 4-6 start at
# position 111 and cover position 150 (offset 39): 1 carries "A" and 2
# carry the alternate "T".
log "Writing SAM reads..."
g80=$(printf 'G%.0s' {1..80})
t80=$(printf 'T%.0s' {1..80})
q80=$(printf 'I%.0s' {1..80})
{
  printf '@HD\tVN:1.6\tSO:unsorted\n'
  printf '@SQ\tSN:chr1\tLN:200\n'
  printf '@RG\tID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\n'
  printf 'read1\t0\tchr1\t1\t60\t80M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n' "$a80" "$q80"
  printf 'read2\t0\tchr1\t1\t60\t80M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n' "$a80" "$q80"
  printf 'read3\t0\tchr1\t1\t60\t80M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n' "$g80" "$q80"
  printf 'read4\t0\tchr1\t111\t60\t80M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n' "$a80" "$q80"
  printf 'read5\t0\tchr1\t111\t60\t80M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n' "$t80" "$q80"
  printf 'read6\t0\tchr1\t111\t60\t80M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n' "$t80" "$q80"
} > "$test_data_dir/reads.sam"
sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.bam"
check_file_exists "$test_data_dir/reads.bam" "sorted BAM file"
check_file_exists "$test_data_dir/reads.bai" "BAM index file"

# The same reads with mapping quality 40, below the tool-default minimum of 50
awk 'BEGIN { OFS = "\t" } /^@/ { print; next } { $5 = 40; print }' "$test_data_dir/reads.sam" > "$test_data_dir/reads_mq40.sam"
sort_and_index_bam "$test_data_dir/reads_mq40.sam" "$test_data_dir/reads_mq40.bam"

# Interval files: BED files for each site and for the whole contig, and a
# VCF with only the site at position 150
printf 'chr1\t49\t50\n' > "$test_data_dir/site_50.bed"
printf 'chr1\t0\t200\n' > "$test_data_dir/whole_contig.bed"
grep -v $'^chr1\t50\t' "$test_data_dir/common_snps.vcf" > "$test_data_dir/site_150.vcf"

# --- Test Case 1: Sites in --variant used as intervals ---
log "Starting TEST 1: Sites in --variant used as intervals"

log "Executing $meta_name with the --variant VCF as --intervals, without indexes or --reference..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups.table"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/pileups.table" "pileup summary table"
check_file_contains "$meta_temp_dir/pileups.table" $'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency' "pileup summary table header"
check_file_contains "$meta_temp_dir/pileups.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts"
check_file_contains "$meta_temp_dir/pileups.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts"
check_file_not_exists "$test_data_dir/common_snps.vcf.idx" "index next to the input VCF (a copy is indexed instead)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: BED intervals and reference ---
log "Starting TEST 2: BED intervals and reference"

log "Executing $meta_name with a BED --intervals file and the reference trio..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/site_50.bed" \
  --reference "$test_data_dir/reference.fasta" \
  --reference_fai "$test_data_dir/reference.fasta.fai" \
  --reference_dict "$test_data_dir/reference.dict" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_bed.table"

log "Validating TEST 2 outputs..."
check_file_contains "$meta_temp_dir/pileups_bed.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts"
check_file_not_contains "$meta_temp_dir/pileups_bed.table" $'chr1\t150\t' "position 150 (not in --intervals)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Pre-built indexes ---
log "Starting TEST 3: Pre-built --variant_index and --intervals_index (.idx and .tbi)"

mkdir -p "$test_data_dir/indexed"
cp "$test_data_dir/common_snps.vcf" "$test_data_dir/indexed/common_snps.vcf"
gatk IndexFeatureFile --input "$test_data_dir/indexed/common_snps.vcf" --verbosity ERROR
check_file_exists "$test_data_dir/indexed/common_snps.vcf.idx" "pre-built Tribble index"
bgzip -c "$test_data_dir/common_snps.vcf" > "$test_data_dir/indexed/common_snps.vcf.gz"
tabix -p vcf "$test_data_dir/indexed/common_snps.vcf.gz"
check_file_exists "$test_data_dir/indexed/common_snps.vcf.gz.tbi" "pre-built tabix index"

log "Executing $meta_name with a .idx --variant_index and a .tbi --intervals_index..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/indexed/common_snps.vcf" \
  --variant_index "$test_data_dir/indexed/common_snps.vcf.idx" \
  --intervals "$test_data_dir/indexed/common_snps.vcf.gz" \
  --intervals_index "$test_data_dir/indexed/common_snps.vcf.gz.tbi" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_indexed.table" > "$meta_temp_dir/test3.log" 2>&1

log "Validating TEST 3 outputs..."
check_file_contains "$meta_temp_dir/pileups_indexed.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts"
check_file_contains "$meta_temp_dir/pileups_indexed.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts"
check_file_not_contains "$meta_temp_dir/test3.log" "IndexFeatureFile" "run log (no indexing when indexes are given)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Interval set rules ---
log "Starting TEST 4: Interval set rules"

log "Executing $meta_name with a BED and a VCF --intervals file (default UNION)..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/site_50.bed" \
  --intervals "$test_data_dir/site_150.vcf" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_union.table"

check_file_contains "$meta_temp_dir/pileups_union.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts (UNION)"
check_file_contains "$meta_temp_dir/pileups_union.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (UNION)"

log "Executing $meta_name with a whole-contig BED and a VCF --intervals file (INTERSECTION)..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/whole_contig.bed" \
  --intervals "$test_data_dir/site_150.vcf" \
  --interval_set_rule INTERSECTION \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_intersection.table"

check_file_contains "$meta_temp_dir/pileups_intersection.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (INTERSECTION)"
check_file_not_contains "$meta_temp_dir/pileups_intersection.table" $'chr1\t50\t' "position 50 (not in the INTERSECTION)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Excluded intervals and allele frequency options ---
log "Starting TEST 5: Excluded intervals and allele frequency options"

log "Executing $meta_name with a BED --exclude_intervals file..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --exclude_intervals "$test_data_dir/site_50.bed" \
  --maximum_population_allele_frequency 0.5 \
  --max_depth_per_sample 100 \
  --read_filter MappingQualityNotZeroReadFilter \
  --output "$meta_temp_dir/pileups_excluded.table"

check_file_contains "$meta_temp_dir/pileups_excluded.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (not excluded)"
check_file_not_contains "$meta_temp_dir/pileups_excluded.table" $'chr1\t50\t' "position 50 (excluded)"

# GATK needs an index for a block-compressed VCF, so this also checks that
# the unindexed .vcf.gz is indexed before the run. The whole contig is used as
# --intervals because GATK fails if -XL removes all of the -L territory.
log "Executing $meta_name with a BED and an unindexed .vcf.gz --exclude_intervals file..."
bgzip -c "$test_data_dir/site_150.vcf" > "$test_data_dir/site_150.vcf.gz"
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/whole_contig.bed" \
  --exclude_intervals "$test_data_dir/site_50.bed" \
  --exclude_intervals "$test_data_dir/site_150.vcf.gz" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_excluded_both.table"

check_file_not_contains "$meta_temp_dir/pileups_excluded_both.table" $'chr1\t50\t' "position 50 (excluded by BED)"
check_file_not_contains "$meta_temp_dir/pileups_excluded_both.table" $'chr1\t150\t' "position 150 (excluded by VCF)"

log "Executing $meta_name with the default --maximum_population_allele_frequency (0.2)..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --output "$meta_temp_dir/pileups_default_af.table"

check_file_not_contains "$meta_temp_dir/pileups_default_af.table" $'chr1\t50\t' "position 50 (AF 0.3 is above the default maximum)"
check_file_not_contains "$meta_temp_dir/pileups_default_af.table" $'chr1\t150\t' "position 150 (AF 0.4 is above the default maximum)"

log "Executing $meta_name with --minimum_population_allele_frequency 0.35..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --minimum_population_allele_frequency 0.35 \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_min_af.table"

check_file_contains "$meta_temp_dir/pileups_min_af.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (AF 0.4)"
check_file_not_contains "$meta_temp_dir/pileups_min_af.table" $'chr1\t50\t' "position 50 (AF 0.3 is below the minimum)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Minimum mapping quality ---
log "Starting TEST 6: Minimum mapping quality"

log "Executing $meta_name on MQ 40 reads with the default --minimum_mapping_quality (50)..."
"$meta_executable" \
  --input "$test_data_dir/reads_mq40.bam" \
  --bai "$test_data_dir/reads_mq40.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_mq_default.table"

# All reads are removed by the tool-default mapping quality read filter, so
# no sites are written
check_file_not_contains "$meta_temp_dir/pileups_mq_default.table" $'chr1\t50\t' "position 50 (MQ 40 reads filtered)"
check_file_not_contains "$meta_temp_dir/pileups_mq_default.table" $'chr1\t150\t' "position 150 (MQ 40 reads filtered)"

log "Executing $meta_name on MQ 40 reads with --minimum_mapping_quality 30..."
"$meta_executable" \
  --input "$test_data_dir/reads_mq40.bam" \
  --bai "$test_data_dir/reads_mq40.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --maximum_population_allele_frequency 0.5 \
  --minimum_mapping_quality 30 \
  --output "$meta_temp_dir/pileups_mq30.table"

check_file_contains "$meta_temp_dir/pileups_mq30.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts (MQ 40 reads kept)"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: --intervals_index count mismatch ---
log "Starting TEST 7: --intervals_index count mismatch fails"

if "$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --bai "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/site_50.bed" \
  --intervals_index "$test_data_dir/indexed/common_snps.vcf.idx" \
  --output "$meta_temp_dir/pileups_bad.table" > "$meta_temp_dir/test7.log" 2>&1; then
  log_error "✗ $meta_name did not fail with an --intervals_index and no VCF --intervals file"
  exit 1
fi
check_file_contains "$meta_temp_dir/test7.log" "must be given once for each VCF file in --intervals" "error message"

log "✅ TEST 7 completed successfully"

print_test_summary "All tests"
