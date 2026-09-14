#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

# --- TEST 1: region query on a bgzipped BED file ---
log "TEST 1: query a region from a bgzipped BED file"
create_test_bed "$test_dir/regions.bed" 3
bgzip -c "$test_dir/regions.bed" >"$test_dir/regions.bed.gz"
tabix -p bed "$test_dir/regions.bed.gz"
check_file_exists "$test_dir/regions.bed.gz.tbi" "BED TBI index"

"$meta_executable" \
  --input "$test_dir/regions.bed.gz" \
  --index "$test_dir/regions.bed.gz.tbi" \
  --region "chr2:1000-2000" \
  --output "$meta_temp_dir/region_result.txt"

check_file_exists "$meta_temp_dir/region_result.txt" "region query result"
check_file_not_empty "$meta_temp_dir/region_result.txt" "region query result"
check_file_contains "$meta_temp_dir/region_result.txt" "region1" "region query result"
log "✅ TEST 1 passed"

# --- TEST 2: list chromosome names from the index ---
log "TEST 2: list chromosome names"
"$meta_executable" \
  --input "$test_dir/regions.bed.gz" \
  --index "$test_dir/regions.bed.gz.tbi" \
  --list_chroms \
  --output "$meta_temp_dir/chroms.txt"

check_file_exists "$meta_temp_dir/chroms.txt" "chromosome list"
check_file_contains "$meta_temp_dir/chroms.txt" "chr2" "chromosome list"
check_file_contains "$meta_temp_dir/chroms.txt" "chr3" "chromosome list"
check_file_contains "$meta_temp_dir/chroms.txt" "chr4" "chromosome list"
log "✅ TEST 2 passed"

# --- TEST 3: region query with header on a bgzipped VCF ---
log "TEST 3: query a region with --print_header on a bgzipped VCF"
cat >"$test_dir/variants.vcf" <<'VCFEOF'
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	50	PASS	.
chr1	500	.	G	C	60	PASS	.
VCFEOF
bgzip -c "$test_dir/variants.vcf" >"$test_dir/variants.vcf.gz"
tabix -p vcf "$test_dir/variants.vcf.gz"
check_file_exists "$test_dir/variants.vcf.gz.tbi" "VCF TBI index"

"$meta_executable" \
  --input "$test_dir/variants.vcf.gz" \
  --index "$test_dir/variants.vcf.gz.tbi" \
  --region "chr1:100-100" \
  --print_header \
  --output "$meta_temp_dir/vcf_result.txt"

check_file_exists "$meta_temp_dir/vcf_result.txt" "VCF region query result"
check_file_contains "$meta_temp_dir/vcf_result.txt" "##fileformat=VCFv4.2" "VCF region query result"
check_file_contains "$meta_temp_dir/vcf_result.txt" "^chr1	100	" "matched VCF record"
log "✅ TEST 3 passed"

# --- TEST 4: only-header mode ---
log "TEST 4: --only_header prints just the header lines"
"$meta_executable" \
  --input "$test_dir/variants.vcf.gz" \
  --index "$test_dir/variants.vcf.gz.tbi" \
  --only_header \
  --output "$meta_temp_dir/header_only.txt"

check_file_exists "$meta_temp_dir/header_only.txt" "header-only result"
check_file_contains "$meta_temp_dir/header_only.txt" "##fileformat=VCFv4.2" "header-only result"
check_file_not_contains "$meta_temp_dir/header_only.txt" "^chr1	100	" "header-only result (no records)"
log "✅ TEST 4 passed"

print_test_summary "$meta_name tests passed"
