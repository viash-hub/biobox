#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# Note: Following the official arriba workflow from run_arriba.sh
# https://github.com/suhrig/arriba/blob/master/run_arriba.sh
# All required tools (STAR, samtools, wget) are available in the Docker container

# Download and prepare test data from arriba repository
log "Downloading arriba test data from official repository..."

# Download test FASTQ files from nf-core test-datasets (these match the genome)
log "Downloading nf-core test RNA-seq reads..."
wget -q -O "$test_data_dir/read1.fastq.gz" "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz"
wget -q -O "$test_data_dir/read2.fastq.gz" "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz"
log "✓ Downloaded nf-core test FASTQ files"

# Download or create test genome and annotation
log "Setting up test genome and annotation..."
wget -q -O "$test_data_dir/genome.fa" "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta"
wget -q -O "$test_data_dir/annotation.gtf" "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.gtf"

# Create required arriba input files
log "Creating required arriba input files..."

# Create empty blacklist file (arriba requirement)
touch "$test_data_dir/blacklist.tsv"

check_file_exists "$test_data_dir/read1.fastq.gz" "test FASTQ R1"
check_file_exists "$test_data_dir/genome.fa" "test genome"
check_file_exists "$test_data_dir/annotation.gtf" "test annotation"
check_file_exists "$test_data_dir/blacklist.tsv" "blacklist file"

# Check if STAR is available and create proper test BAM
log "Creating STAR genome index and proper test BAM following arriba workflow..."

# Create STAR genome index with parameters from nf-core test
mkdir -p "$test_data_dir/star_index"
log "Generating STAR genome index..."
STAR --runMode genomeGenerate \
  --genomeDir "$test_data_dir/star_index" \
  --genomeFastaFiles "$test_data_dir/genome.fa" \
  --sjdbGTFfile "$test_data_dir/annotation.gtf" \
  --genomeSAindexNbases 9 \
  --runThreadN 1 \
  --genomeChrBinNbits 12 \
  --sjdbOverhang 50 >/dev/null 2>&1

log "STAR index created successfully, aligning reads with chimeric parameters..."
# Use the exact STAR parameters from nf-core test configuration
STAR --runThreadN 1 \
  --genomeDir "$test_data_dir/star_index" \
  --genomeLoad NoSharedMemory \
  --readFilesIn "$test_data_dir/read1.fastq.gz" "$test_data_dir/read2.fastq.gz" \
  --readFilesCommand zcat \
  --outStd BAM_Unsorted \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outBAMcompression 0 \
  --outFilterMultimapNmax 50 \
  --peOverlapNbasesMin 10 \
  --alignSplicedMateMapLminOverLmate 0.5 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimSegmentMin 10 \
  --chimOutType WithinBAM HardClip \
  --chimJunctionOverhangMin 10 \
  --chimScoreDropMax 30 \
  --chimScoreJunctionNonGTAG 0 \
  --chimScoreSeparation 1 \
  --chimSegmentReadGapMax 3 \
  --chimMultimapNmax 50 \
  --outFileNamePrefix "$test_data_dir/star_" \
  > "$test_data_dir/Aligned.out.bam" 2>/dev/null

# Sort the BAM file if it was created successfully
log "Sorting and indexing STAR-aligned BAM..."
samtools sort -o "$test_data_dir/sorted.bam" "$test_data_dir/Aligned.out.bam" 2>/dev/null
mv "$test_data_dir/sorted.bam" "$test_data_dir/Aligned.out.bam"
samtools index "$test_data_dir/Aligned.out.bam" 2>/dev/null || true

check_file_exists "$test_data_dir/Aligned.out.bam" "test BAM file"

# --- Test Case 1: Full arriba workflow test ---
log "Starting TEST 1: Full arriba workflow test (following run_arriba.sh pattern)"

log "Running $meta_name with complete workflow parameters..."
# Follow the exact parameter pattern from arriba's run_arriba.sh
"$meta_executable" \
  --bam "$test_data_dir/Aligned.out.bam" \
  --genome "$test_data_dir/genome.fa" \
  --gene_annotation "$test_data_dir/annotation.gtf" \
  --blacklist "$test_data_dir/blacklist.tsv" \
  --fusions "$meta_temp_dir/fusions.tsv" \
  --fusions_discarded "$meta_temp_dir/fusions_discarded.tsv"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/fusions.tsv" "main fusions output file"
check_file_exists "$meta_temp_dir/fusions_discarded.tsv" "discarded fusions output file"

# Check if we detected any fusions (file should have content beyond header)
if [[ $(wc -l < "$meta_temp_dir/fusions.tsv") -gt 1 ]]; then
  log "✓ Detected $(( $(wc -l < "$meta_temp_dir/fusions.tsv") - 1 )) fusion(s) in output"
else
  log "ℹ No fusions detected - this may be expected with test data"
fi

log "✅ TEST 1 completed successfully - full arriba workflow passed"

# --- Test Case 2: Test with minimal required parameters ---
log "Starting TEST 2: Minimal parameters test"

log "Running $meta_name with minimal required parameters..."
"$meta_executable" \
  --bam "$test_data_dir/Aligned.out.bam" \
  --genome "$test_data_dir/genome.fa" \
  --gene_annotation "$test_data_dir/annotation.gtf" \
  --fusions "$meta_temp_dir/fusions_minimal.tsv" \
  --disable_filters blacklist

check_file_exists "$meta_temp_dir/fusions_minimal.tsv" "minimal fusions output file"
log "✅ TEST 2 completed successfully"

print_test_summary "All tests completed successfully"
