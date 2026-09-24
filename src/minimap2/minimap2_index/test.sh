#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

cat > "$test_dir/ref.fasta" <<'EOF'
>chr1
ACTGATCGATCGATCGATCGATCGATCGATCGATCGATCGACTATCGATCGATCGATCGA
EOF
check_file_exists "$test_dir/ref.fasta" "test reference FASTA"

# --- TEST 1: create a minimap2 index ---
log "TEST 1: create a minimap2 index from a FASTA"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --output "$meta_temp_dir/ref.mmi"

check_file_exists "$meta_temp_dir/ref.mmi" "minimap2 index"
check_file_not_empty "$meta_temp_dir/ref.mmi" "minimap2 index"

# minimap2 index files start with the magic bytes "MMI\2"
if ! head -c 3 "$meta_temp_dir/ref.mmi" | grep -q "MMI"; then
  log_error "output does not look like a minimap2 index (missing MMI magic bytes)"
  exit 1
fi
log "- minimap2 index magic bytes present"

log "OK: TEST 1 passed"

# --- TEST 2: --preset selects the indexing parameters ---
# map-hifi implies k=19, w=19; the default is k=15, w=10. The parameters are
# baked into the .mmi, so this is what makes a pre-built index usable for a
# preset other than the long-read default.
log "TEST 2: --preset map-hifi bakes k=19/w=19 into the index"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --preset map-hifi \
  --output "$meta_temp_dir/ref_hifi.mmi" \
  > "$meta_temp_dir/hifi_index.log" 2>&1

check_file_exists "$meta_temp_dir/ref_hifi.mmi" "map-hifi index"
check_file_not_empty "$meta_temp_dir/ref_hifi.mmi" "map-hifi index"
check_file_matches_regex "$meta_temp_dir/hifi_index.log" \
  "kmer size: 19; skip: 19" "map-hifi indexing parameters (k=19, w=19)"
log "OK: TEST 2 passed"

# --- TEST 3: explicit --kmer_size / --window_size override the preset ---
log "TEST 3: --kmer_size and --window_size override the preset"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --preset map-hifi \
  --kmer_size 21 \
  --window_size 11 \
  --output "$meta_temp_dir/ref_k21.mmi" \
  > "$meta_temp_dir/k21_index.log" 2>&1

check_file_not_empty "$meta_temp_dir/ref_k21.mmi" "k=21 index"
check_file_matches_regex "$meta_temp_dir/k21_index.log" \
  "kmer size: 21; skip: 11" "overridden indexing parameters (k=21, w=11)"
log "OK: TEST 3 passed"

print_test_summary "$meta_name tests passed"
