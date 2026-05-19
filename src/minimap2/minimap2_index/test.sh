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
log "✓ minimap2 index magic bytes present"

log "✅ TEST 1 passed"

print_test_summary "$meta_name tests passed"
