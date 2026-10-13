#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Set up the test environment and the test data directory
gatk4_test_init

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
reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

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


# The same reads with other read names in another read group of the same
# sample, as for a second lane
awk 'BEGIN { OFS = "\t" }
  /^@RG/ { sub(/ID:rg1/, "ID:rg2"); print; next }
  /^@/ { print; next }
  { $1 = "lane2_" $1; $12 = "RG:Z:rg2"; print }' "$test_data_dir/reads.sam" > "$test_data_dir/reads_rg2.sam"
sort_and_index_bam "$test_data_dir/reads_rg2.sam" "$test_data_dir/reads_rg2.bam"
check_file_exists "$test_data_dir/reads_rg2.bai" "second BAM index"

# The reads as a CRAM with its index in another directory
create_test_cram "$test_data_dir/reads.bam" "$test_data_dir/reference.fasta" \
  "$test_data_dir/cram/reads.cram" "$test_data_dir/cram_index/reads_index.crai"

# Copies of the common SNPs VCF with pre-built indexes (.idx and .tbi)
mkdir -p "$test_data_dir/indexed"
cp "$test_data_dir/common_snps.vcf" "$test_data_dir/indexed/common_snps.vcf"
gatk IndexFeatureFile --input "$test_data_dir/indexed/common_snps.vcf" --verbosity ERROR
check_file_exists "$test_data_dir/indexed/common_snps.vcf.idx" "pre-built Tribble index"
bgzip -c "$test_data_dir/common_snps.vcf" > "$test_data_dir/indexed/common_snps.vcf.gz"
tabix -p vcf "$test_data_dir/indexed/common_snps.vcf.gz"
check_file_exists "$test_data_dir/indexed/common_snps.vcf.gz.tbi" "pre-built tabix index"

# Unindexed bgzipped copies of the position 150 VCF with the .vcf.gz and
# .vcf.bgz extensions, and a .vcf.gz compressed with gzip (not bgzip), which
# cannot be indexed
bgzip -c "$test_data_dir/site_150.vcf" > "$test_data_dir/site_150.vcf.gz"
bgzip -c "$test_data_dir/site_150.vcf" > "$test_data_dir/site_150.vcf.bgz"
gzip -c "$test_data_dir/common_snps.vcf" > "$test_data_dir/plain_gzip.vcf.gz"

# Positions 146-148, which end 2bp before the site at position 150
printf 'chr1\t145\t148\n' > "$test_data_dir/near_150.bed"

# --- Test Case 1: Basic run ---
log "Starting TEST 1: Basic run (sites in --variant used as intervals)"

log "Executing $meta_name with the --variant VCF as --intervals, without indexes or --reference..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
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

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: CRAM and BAM inputs with the reference, and pre-built indexes"

# The doubled counts prove that the reads of both inputs (the CRAM and the
# second BAM) are counted
log "Executing $meta_name with a CRAM (--read_index in another directory) and a BAM --input, and the reference trio..."
"$meta_executable" \
  --input "$test_data_dir/cram/reads.cram" \
  --read_index "$test_data_dir/cram_index/reads_index.crai" \
  --input "$test_data_dir/reads_rg2.bam" \
  --read_index "$test_data_dir/reads_rg2.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  "${reference_args[@]}" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_two_inputs.table"

check_file_contains "$meta_temp_dir/pileups_two_inputs.table" $'chr1\t50\t4\t2\t0\t0.3' "position 50 counts (reads of both inputs)"
check_file_contains "$meta_temp_dir/pileups_two_inputs.table" $'chr1\t150\t2\t4\t0\t0.4' "position 150 counts (reads of both inputs)"

log "Executing $meta_name with a .idx --variant_index and a .tbi --intervals_index..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/indexed/common_snps.vcf" \
  --variant_index "$test_data_dir/indexed/common_snps.vcf.idx" \
  --intervals "$test_data_dir/indexed/common_snps.vcf.gz" \
  --intervals_index "$test_data_dir/indexed/common_snps.vcf.gz.tbi" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_indexed.table" > "$meta_temp_dir/test2_indexed.log" 2>&1

check_file_contains "$meta_temp_dir/pileups_indexed.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts"
check_file_contains "$meta_temp_dir/pileups_indexed.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts"
check_file_not_contains "$meta_temp_dir/test2_indexed.log" "IndexFeatureFile" "run log (no indexing when indexes are given)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Interval set rules, excluded intervals and interval padding"

log "Executing $meta_name with a BED and an unindexed .vcf.bgz --intervals file (default UNION)..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/site_50.bed" \
  --intervals "$test_data_dir/site_150.vcf.bgz" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_union.table" 2>&1 | tee "$meta_temp_dir/test3_union.log"

check_file_contains "$meta_temp_dir/pileups_union.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts (UNION, BED)"
check_file_contains "$meta_temp_dir/pileups_union.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (UNION, .vcf.bgz)"
check_file_contains "$meta_temp_dir/test3_union.log" "Warning: no index was provided for '.*site_150.vcf.bgz'" "run log (.vcf.bgz staged as a VCF)"

log "Executing $meta_name with a whole-contig BED and a VCF --intervals file (INTERSECTION)..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/whole_contig.bed" \
  --intervals "$test_data_dir/site_150.vcf" \
  --interval_set_rule INTERSECTION \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_intersection.table"

check_file_contains "$meta_temp_dir/pileups_intersection.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (INTERSECTION)"
check_file_not_contains "$meta_temp_dir/pileups_intersection.table" $'chr1\t50\t' "position 50 (not in the INTERSECTION)"

# GATK needs an index for a block-compressed VCF, so this also checks that
# the unindexed .vcf.gz is indexed before the run. The whole contig is used as
# --intervals because GATK fails if -XL removes all of the -L territory.
log "Executing $meta_name with an unindexed .vcf.gz --exclude_intervals file, --max_depth_per_sample and --read_filter..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/whole_contig.bed" \
  --exclude_intervals "$test_data_dir/site_150.vcf.gz" \
  --maximum_population_allele_frequency 0.5 \
  --max_depth_per_sample 100 \
  --read_filter MappingQualityNotZeroReadFilter \
  --output "$meta_temp_dir/pileups_excluded.table" 2>&1 | tee "$meta_temp_dir/test3_excluded.log"

check_file_contains "$meta_temp_dir/pileups_excluded.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts (not excluded)"
check_file_not_contains "$meta_temp_dir/pileups_excluded.table" $'chr1\t150\t' "position 150 (excluded by the .vcf.gz)"
check_file_contains "$meta_temp_dir/test3_excluded.log" "max-depth-per-sample 100" "run log command line (--max_depth_per_sample)"
check_file_contains "$meta_temp_dir/test3_excluded.log" "read-filter MappingQualityNotZeroReadFilter" "run log command line (--read_filter)"

log "Executing $meta_name with an --intervals file next to position 150, without padding..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/near_150.bed" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_no_padding.table"

check_file_not_contains "$meta_temp_dir/pileups_no_padding.table" $'chr1\t150\t' "position 150 (outside the interval)"

log "Executing $meta_name with the same --intervals file and --interval_padding 2..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/near_150.bed" \
  --interval_padding 2 \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_padding.table"

check_file_contains "$meta_temp_dir/pileups_padding.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (in the padded interval)"

log "Executing $meta_name with the same file as a BED --exclude_intervals file and --interval_exclusion_padding 2..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/whole_contig.bed" \
  --exclude_intervals "$test_data_dir/near_150.bed" \
  --interval_exclusion_padding 2 \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_exclusion_padding.table"

check_file_contains "$meta_temp_dir/pileups_exclusion_padding.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts (not excluded)"
check_file_not_contains "$meta_temp_dir/pileups_exclusion_padding.table" $'chr1\t150\t' "position 150 (in the padded excluded interval)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Tool options ---
log "Starting TEST 4: Allele frequency and mapping quality options"

log "Executing $meta_name with the default --maximum_population_allele_frequency (0.2)..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --output "$meta_temp_dir/pileups_default_af.table"

check_file_not_contains "$meta_temp_dir/pileups_default_af.table" $'chr1\t50\t' "position 50 (AF 0.3 is above the default maximum)"
check_file_not_contains "$meta_temp_dir/pileups_default_af.table" $'chr1\t150\t' "position 150 (AF 0.4 is above the default maximum)"

log "Executing $meta_name with --minimum_population_allele_frequency 0.35..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --minimum_population_allele_frequency 0.35 \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_min_af.table"

check_file_contains "$meta_temp_dir/pileups_min_af.table" $'chr1\t150\t1\t2\t0\t0.4' "position 150 counts (AF 0.4)"
check_file_not_contains "$meta_temp_dir/pileups_min_af.table" $'chr1\t50\t' "position 50 (AF 0.3 is below the minimum)"

log "Executing $meta_name on MQ 40 reads with the default --minimum_mapping_quality (50)..."
"$meta_executable" \
  --input "$test_data_dir/reads_mq40.bam" \
  --read_index "$test_data_dir/reads_mq40.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --maximum_population_allele_frequency 0.5 \
  --output "$meta_temp_dir/pileups_mq_default.table"

# All reads are removed by the tool-default mapping quality read filter, so
# no sites are written
check_file_not_contains "$meta_temp_dir/pileups_mq_default.table" $'chr1\t50\t' "position 50 (MQ 40 reads filtered)"
check_file_not_contains "$meta_temp_dir/pileups_mq_default.table" $'chr1\t150\t' "position 150 (MQ 40 reads filtered)"

log "Executing $meta_name on MQ 60 reads with --maximum_mapping_quality 50..."
"$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --maximum_population_allele_frequency 0.5 \
  --minimum_mapping_quality 30 \
  --maximum_mapping_quality 50 \
  --output "$meta_temp_dir/pileups_max_mq.table"

# All MQ 60 reads are removed, so no sites are written
check_file_not_contains "$meta_temp_dir/pileups_max_mq.table" $'chr1\t50\t' "position 50 (MQ 60 reads filtered)"
check_file_not_contains "$meta_temp_dir/pileups_max_mq.table" $'chr1\t150\t' "position 150 (MQ 60 reads filtered)"

# Only --interval_strings (no --intervals file) is given, and it does not
# include position 150
log "Executing $meta_name on MQ 40 reads with --minimum_mapping_quality 30, --maximum_mapping_quality 50 and only --interval_strings..."
"$meta_executable" \
  --input "$test_data_dir/reads_mq40.bam" \
  --read_index "$test_data_dir/reads_mq40.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --interval_strings "chr1:1-100" \
  --maximum_population_allele_frequency 0.5 \
  --minimum_mapping_quality 30 \
  --maximum_mapping_quality 50 \
  --output "$meta_temp_dir/pileups_mq40_strings.table"

check_file_contains "$meta_temp_dir/pileups_mq40_strings.table" $'chr1\t50\t2\t1\t0\t0.3' "position 50 counts (MQ 40 reads kept, in the interval string)"
check_file_not_contains "$meta_temp_dir/pileups_mq40_strings.table" $'chr1\t150\t' "position 150 (outside the interval string)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Argument errors ---
log "Starting TEST 5: Argument errors"

check_fails_with "must be given once for each VCF file in --intervals" "$meta_temp_dir/test5_intervals_index.log" \
  "$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/site_50.bed" \
  --intervals_index "$test_data_dir/indexed/common_snps.vcf.idx" \
  --output "$meta_temp_dir/pileups_bad.table"

check_fails_with "Error: could not index '" "$meta_temp_dir/test5_gzip.log" \
  "$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/plain_gzip.vcf.gz" \
  --intervals "$test_data_dir/site_50.bed" \
  --output "$meta_temp_dir/pileups_gzip.table"
check_file_contains "$meta_temp_dir/test5_gzip.log" "Error: could not index '.*plain_gzip.vcf.gz'" "error message"
check_file_not_contains "$meta_temp_dir/test5_gzip.log" "GetPileupSummaries done" "run log (GATK did not run)"

check_fails_with "Error: --intervals or --interval_strings is required." "$meta_temp_dir/test5_no_intervals.log" \
  "$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --output "$meta_temp_dir/pileups_no_intervals.table"
check_file_not_contains "$meta_temp_dir/test5_no_intervals.log" "GetPileupSummaries done" "run log (GATK did not run)"

check_fails_with "Error: --reference is required for CRAM input or output." "$meta_temp_dir/test5_cram_no_ref.log" \
  "$meta_executable" \
  --input "$test_data_dir/cram/reads.cram" \
  --read_index "$test_data_dir/cram_index/reads_index.crai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --output "$meta_temp_dir/pileups_cram_no_ref.table"
check_file_not_contains "$meta_temp_dir/test5_cram_no_ref.log" "GetPileupSummaries done" "run log (GATK did not run)"

check_fails_with "Error: --read_index must be given once for each file in --input" "$meta_temp_dir/test5_read_index.log" \
  "$meta_executable" \
  --input "$test_data_dir/reads.bam" \
  --input "$test_data_dir/reads_rg2.bam" \
  --read_index "$test_data_dir/reads.bai" \
  --variant "$test_data_dir/common_snps.vcf" \
  --intervals "$test_data_dir/common_snps.vcf" \
  --output "$meta_temp_dir/pileups_bad_index.table"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests"
