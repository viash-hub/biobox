#!/bin/bash

## VIASH START
## VIASH END

# source the helpers & setup test env
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Test for $meta_name starting ..."

# create and validate test-data
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir" "$meta_temp_dir/index"
create_test_fasta "$test_data_dir/22_20-21M.fa" 1 10000
check_file_exists "$test_data_dir/22_20-21M.fa" "test reference FASTA"

# --- TEST_1: create a hisat2 index ---
log "Starting TEST_1: Create a hisat2 index from a reference FASTA"
"$meta_executable" \
  --reference "$test_data_dir/22_20-21M.fa" \
  --index_dir "$meta_temp_dir/index" \
  --index_prefix genome

for ext in 1.ht2 2.ht2 3.ht2 4.ht2 5.ht2 6.ht2 7.ht2 8.ht2; do
  check_file_exists "$meta_temp_dir/index/genome.$ext" "index file genome.$ext"
done

log "✓ TEST_1 completed successfully: 8 index files were created"

# --- TEST 2: splice-aware index with --ss and --exon ---
log "Starting TEST_2: building a splice-aware (graph) index with --ss and --exon"

# Splice sites: chrom, donor (0-based), acceptor (0-based), strand
printf "seq1\t1000\t2000\t+\n" > "$test_data_dir/splicesites.txt"
# Exons: chrom, start (0-based, inclusive), end (0-based, exclusive)
printf "seq1\t0\t1000\nseq1\t2000\t3000\n" > "$test_data_dir/exons.txt"

build_log="$meta_temp_dir/build_splice.log"
"$meta_executable" \
  --reference "$test_data_dir/22_20-21M.fa" \
  --index_dir "$meta_temp_dir/index_splice" \
  --index_prefix genome_splice \
  --ss "$test_data_dir/splicesites.txt" \
  --exon "$test_data_dir/exons.txt" \
  2>&1 | tee "$build_log"

for ext in 1.ht2 2.ht2 3.ht2 4.ht2 5.ht2 6.ht2 7.ht2 8.ht2; do
  check_file_exists "$meta_temp_dir/index_splice/genome_splice.$ext" "splice-aware index file genome_splice.$ext"
done

# check build log on "LinearFM: No"
check_file_matches_regex "$build_log" "linearFM: No" "graph FM-index indicator"

log "✓ TEST_2 completed successfully: splice-aware graph index created"

print_test_summary "$meta_name tests"