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
>seq1
ACTGATCGATCGATCGATCGATCGATCGATCGATCGATCGACTATCGATCGATCGATCGA
EOF

cat > "$test_dir/reads.fastq" <<'EOF'
@read1
ACTGATCGATCGATCGATCGATCAAAAGATCGATCGATCGACTATCGATCGATCGATCGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

check_file_exists "$test_dir/ref.fasta" "reference FASTA"
check_file_exists "$test_dir/reads.fastq" "query FASTQ"

# --- TEST 1: PAF output ---
log "TEST 1: alignment to PAF"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --output "$meta_temp_dir/out.paf"

check_file_exists "$meta_temp_dir/out.paf" "PAF output"
check_file_not_empty "$meta_temp_dir/out.paf" "PAF output"
check_file_contains "$meta_temp_dir/out.paf" "seq1" "PAF alignment to seq1"
log "✅ TEST 1 passed"

# --- TEST 2: sorted + indexed BAM output ---
log "TEST 2: alignment to sorted BAM"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/out.bam"

check_file_exists "$meta_temp_dir/out.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/out.bam" "BAM output"
check_file_exists "$meta_temp_dir/out.bam.bai" "BAM index"

if ! samtools view "$meta_temp_dir/out.bam" | grep -q "seq1"; then
  log_error "BAM alignment to seq1 missing"
  exit 1
fi
log "✓ BAM contains alignment to seq1"
log "✅ TEST 2 passed"

print_test_summary "$meta_name tests passed"
