#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_igv"

# Create test data following documentation guidelines
log "Creating test data..."

# Create basic intervals file with name field
cat > "$meta_temp_dir/intervals.bed" << 'EOF'
chr1	1000	2000	region1	100	+
chr1	5000	6000	region2	200	-
chr2	10000	11000	region3	150	+
chr2	20000	21000	region4	300	-
chr3	30000	31000	region5	250	+
EOF

# Create intervals without name field
cat > "$meta_temp_dir/simple.bed" << 'EOF'
chr1	2000	3000
chr1	7000	8000
chr2	15000	16000
EOF

# Create GFF test file
cat > "$meta_temp_dir/features.gff" << 'EOF'
##gff-version 3
chr1	source	gene	1500	2500	.	+	.	ID=gene1;Name=TestGene1
chr1	source	exon	1500	1800	.	+	.	ID=exon1;Parent=gene1
chr1	source	exon	2200	2500	.	+	.	ID=exon2;Parent=gene1
chr2	source	gene	12000	13000	.	-	.	ID=gene2;Name=TestGene2
EOF

# Create mock IGV session file
cat > "$meta_temp_dir/session.xml" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="hg19" locus="chr1:1000-2000">
  <Files>
      <DataFile name="Test Track" path="/path/to/test.bam"/>
  </Files>
</Session>
EOF

# TEST 1: Basic IGV batch script generation
log "Starting TEST 1: Basic IGV batch script generation"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --output "$meta_temp_dir/basic_script.txt"

check_file_exists "$meta_temp_dir/basic_script.txt" "basic IGV script"
check_file_not_empty "$meta_temp_dir/basic_script.txt" "basic IGV script"

# Check that script contains expected IGV commands
if grep -q "snapshot" "$meta_temp_dir/basic_script.txt"; then
  log "✓ basic script contains snapshot commands: $meta_temp_dir/basic_script.txt"
else
  log "✗ basic script missing snapshot commands: $meta_temp_dir/basic_script.txt"
  exit 1
fi

# Check that script contains goto commands for each region
region_count=$(grep -c "goto" "$meta_temp_dir/basic_script.txt" || true)
if [ "$region_count" -eq 5 ]; then
  log "✓ basic script contains expected number of goto commands (5): $meta_temp_dir/basic_script.txt"
else
  log "✗ basic script has unexpected goto command count ($region_count, expected 5): $meta_temp_dir/basic_script.txt"
  exit 1
fi

log "✅ TEST 1 completed successfully"

# TEST 2: IGV script with output path and image format
log "Starting TEST 2: IGV script with custom output path and format"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --output_path "/custom/path/images/" \
  --image_format "svg" \
  --output "$meta_temp_dir/custom_script.txt"

check_file_exists "$meta_temp_dir/custom_script.txt" "custom IGV script"
check_file_not_empty "$meta_temp_dir/custom_script.txt" "custom IGV script"

# Check for custom output path in script
if grep -q "/custom/path/images/" "$meta_temp_dir/custom_script.txt"; then
  log "✓ custom script contains specified output path: $meta_temp_dir/custom_script.txt"
else
  log "✗ custom script missing specified output path: $meta_temp_dir/custom_script.txt"
  exit 1
fi

# Check for SVG format specification
if grep -q "svg" "$meta_temp_dir/custom_script.txt"; then
  log "✓ custom script specifies SVG format: $meta_temp_dir/custom_script.txt"
else
  log "✗ custom script missing SVG format: $meta_temp_dir/custom_script.txt"
  exit 1
fi

log "✅ TEST 2 completed successfully"

# TEST 3: IGV script with session file loading
log "Starting TEST 3: IGV script with session file"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --session_file "$meta_temp_dir/session.xml" \
  --output "$meta_temp_dir/session_script.txt"

check_file_exists "$meta_temp_dir/session_script.txt" "session IGV script"
check_file_not_empty "$meta_temp_dir/session_script.txt" "session IGV script"

# Check for session loading command
if grep -q "session.xml" "$meta_temp_dir/session_script.txt"; then
  log "✓ session script contains session file reference: $meta_temp_dir/session_script.txt"
else
  log "✗ session script missing session file reference: $meta_temp_dir/session_script.txt"
  exit 1
fi

log "✅ TEST 3 completed successfully"

# TEST 4: IGV script with read sorting and collapse
log "Starting TEST 4: IGV script with read display options"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --sort_reads "position" \
  --collapse_reads \
  --output "$meta_temp_dir/display_script.txt"

check_file_exists "$meta_temp_dir/display_script.txt" "display options IGV script"
check_file_not_empty "$meta_temp_dir/display_script.txt" "display options IGV script"

# Check for sorting command
if grep -q "sort" "$meta_temp_dir/display_script.txt"; then
  log "✓ display script contains sorting commands: $meta_temp_dir/display_script.txt"
else
  log "✗ display script missing sorting commands: $meta_temp_dir/display_script.txt"
  exit 1
fi

log "✅ TEST 4 completed successfully"

# TEST 5: IGV script with flanking regions and name-based filenames
log "Starting TEST 5: IGV script with flanking and named files"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --flank_size 500 \
  --use_name \
  --output "$meta_temp_dir/flanked_script.txt"

check_file_exists "$meta_temp_dir/flanked_script.txt" "flanked IGV script"
check_file_not_empty "$meta_temp_dir/flanked_script.txt" "flanked IGV script"

# Check for expanded regions (should include flanking) - chr1:5000-6000 with 500bp flanking = chr1:4500-6500
if grep -q "4500-6500" "$meta_temp_dir/flanked_script.txt"; then
  log "✓ flanked script contains expanded regions: $meta_temp_dir/flanked_script.txt"
else
  log "✗ flanked script missing expanded regions: $meta_temp_dir/flanked_script.txt"
  cat "$meta_temp_dir/flanked_script.txt" >&2
  exit 1
fi

log "✅ TEST 5 completed successfully"

# TEST 6: IGV script with GFF input
log "Starting TEST 6: IGV script with GFF input"
"$meta_executable" \
  --input "$meta_temp_dir/features.gff" \
  --output "$meta_temp_dir/gff_script.txt"

check_file_exists "$meta_temp_dir/gff_script.txt" "GFF IGV script"
check_file_not_empty "$meta_temp_dir/gff_script.txt" "GFF IGV script"

# Should contain regions from GFF file
gene_count=$(grep -c "goto" "$meta_temp_dir/gff_script.txt" || true)
if [ "$gene_count" -ge 2 ]; then
  log "✓ GFF script contains expected regions (≥2): $meta_temp_dir/gff_script.txt"
else
  log "✗ GFF script has too few regions ($gene_count, expected ≥2): $meta_temp_dir/gff_script.txt"
  exit 1
fi

log "✅ TEST 6 completed successfully"

# TEST 7: IGV script with minimal BED input (no name field)
log "Starting TEST 7: IGV script with simple BED input"
"$meta_executable" \
  --input "$meta_temp_dir/simple.bed" \
  --image_format "png" \
  --output "$meta_temp_dir/simple_script.txt"

check_file_exists "$meta_temp_dir/simple_script.txt" "simple BED IGV script"
check_file_not_empty "$meta_temp_dir/simple_script.txt" "simple BED IGV script"

# Should work with 3-column BED format
simple_count=$(grep -c "goto" "$meta_temp_dir/simple_script.txt" || true)
if [ "$simple_count" -eq 3 ]; then
  log "✓ simple script handles 3-column BED correctly (3 regions): $meta_temp_dir/simple_script.txt"
else
  log "✗ simple script region count mismatch ($simple_count, expected 3): $meta_temp_dir/simple_script.txt"
  exit 1
fi

log "✅ TEST 7 completed successfully"

log "All tests completed successfully!"
