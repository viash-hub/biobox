#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

# --- Test data -------------------------------------------------------------
# A 20 kb reference of deterministic pseudo-random sequence (MINSTD LCG in awk;
# the engine image has no python3). A short, low-complexity toy sequence is not
# enough: with the long-read presets (k=19, w=19) it yields no minimizers, so
# an index "built" from it proves nothing.
awk 'BEGIN {
  bases = "ACGT"
  x = 42
  line = ""
  print ">seq1"
  for (i = 0; i < 20000; i++) {
    x = (x * 16807) % 2147483647
    line = line substr(bases, int(x / 2147483647 * 4) + 1, 1)
    if (length(line) == 70) { print line; line = "" }
  }
  if (length(line) > 0) print line
}' > "$test_dir/ref.fasta"
check_file_not_empty "$test_dir/ref.fasta" "test reference FASTA"

# --- Helpers ---------------------------------------------------------------
# A .mmi starts with the magic "MMI\2" followed by the uint32 fields
# w, k, b, n_seq, flags. Reading them from the index itself proves which
# parameters were baked in, rather than trusting minimap2's log.
# Prints "w k n_seq hpc", where hpc is flag bit 0 (MM_I_HPC).
mmi_params() {
  # od wraps after four values, so join its lines before picking fields
  od -An -tu4 -j4 -N20 "$1" | tr -s ' \n' '  ' | awk '{ printf "%s %s %s %s\n", $1, $2, $4, $5 % 2 }'
}

check_mmi() {
  local mmi="$1" expected="$2" label="$3" actual
  check_file_not_empty "$mmi" "$label"
  if [[ "$(head -c 3 "$mmi")" != "MMI" ]]; then
    log_error "$label does not look like a minimap2 index (missing MMI magic bytes)"
    exit 1
  fi
  actual=$(mmi_params "$mmi")
  if [[ "$actual" != "$expected" ]]; then
    log_error "$label: index header 'w k n_seq hpc' is '$actual', expected '$expected'"
    exit 1
  fi
  log "- $label: index header 'w k n_seq hpc' = $actual"
}

# A negative test that only checks the exit code passes on any failure at all,
# so each one is pinned to the error message.
check_rejected() {
  local log_file="$1" message="$2" label="$3"
  shift 3
  if "$meta_executable" "$@" > "$log_file" 2>&1; then
    log_error "$label should have failed"
    exit 1
  fi
  check_file_contains "$log_file" "$message" "error message ($label)"
}

# --- TEST 1: default index ---
log "TEST 1: create a minimap2 index from a FASTA"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --output "$meta_temp_dir/ref.mmi"

# minimap2 defaults: w=10, k=15, one sequence, no homopolymer compression
check_mmi "$meta_temp_dir/ref.mmi" "10 15 1 0" "default index"
log "OK: TEST 1 passed"

# --- TEST 2: --preset selects the indexing parameters ---
# map-hifi implies k=19, w=19. The parameters are baked into the .mmi, so this
# is what makes a pre-built index usable for a preset other than the default.
log "TEST 2: --preset map-hifi bakes k=19/w=19 into the index"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --preset map-hifi \
  --output "$meta_temp_dir/ref_hifi.mmi"

check_mmi "$meta_temp_dir/ref_hifi.mmi" "19 19 1 0" "map-hifi index"
log "OK: TEST 2 passed"

# --- TEST 3: explicit --kmer_size / --window_size override the preset ---
log "TEST 3: --kmer_size and --window_size override the preset"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --preset map-hifi \
  --kmer_size 21 \
  --window_size 11 \
  --output "$meta_temp_dir/ref_k21.mmi"

check_mmi "$meta_temp_dir/ref_k21.mmi" "11 21 1 0" "k=21/w=11 index"
log "OK: TEST 3 passed"

# --- TEST 4: --homopolymer_compressed sets the HPC flag ---
log "TEST 4: --homopolymer_compressed"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --homopolymer_compressed \
  --output "$meta_temp_dir/ref_hpc.mmi"

check_mmi "$meta_temp_dir/ref_hpc.mmi" "10 15 1 1" "homopolymer-compressed index"
log "OK: TEST 4 passed"

# --- TEST 5: index under a cpu allocation ---
# ---cpus is what populates meta_cpus, the only way the -t branch runs at all.
# Threading must not change the index.
log "TEST 5: index under ---cpus 2"
"$meta_executable" \
  --input "$test_dir/ref.fasta" \
  --output "$meta_temp_dir/ref_threads.mmi" \
  ---cpus 2 \
  > "$meta_temp_dir/ref_threads.log" 2>&1

# An identical index would also result if -t were never passed on, so check
# the command line minimap2 echoes.
check_file_contains "$meta_temp_dir/ref_threads.log" "CMD: minimap2 -t 2" "index log (minimap2 -t)"
if ! cmp -s "$meta_temp_dir/ref.mmi" "$meta_temp_dir/ref_threads.mmi"; then
  log_error "Index built with ---cpus 2 differs from the single-threaded index"
  exit 1
fi
log "- byte-identical to the default index"
log "OK: TEST 5 passed"

# --- TEST 6: out-of-range k / w are rejected ---
# minimap2 itself does not reject k > 28 cleanly: it aborts on an assertion in
# mm_sketch (exit 134). The config bounds catch it before minimap2 runs.
log "TEST 6: out-of-range --kmer_size / --window_size are rejected"
check_rejected "$meta_temp_dir/bad_k.log" "kmer_size" "kmer_size 30" \
  --input "$test_dir/ref.fasta" --kmer_size 30 --output "$meta_temp_dir/bad_k.mmi"
check_file_not_contains "$meta_temp_dir/bad_k.log" "Assertion" "kmer_size 30 log (minimap2 assertion)"
check_rejected "$meta_temp_dir/bad_w.log" "window_size" "window_size 256" \
  --input "$test_dir/ref.fasta" --window_size 256 --output "$meta_temp_dir/bad_w.mmi"
log "OK: TEST 6 passed"

print_test_summary "$meta_name tests passed"
