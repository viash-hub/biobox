#!/bin/bash

set -e

## VIASH START
meta_executable="target/docker/bd_rhapsody/bd_rhapsody_sequence_analysis/bd_rhapsody_sequence_analysis"
meta_resources_dir="src/bd_rhapsody/bd_rhapsody_sequence_analysis"
## VIASH END

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_empty() {
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
assert_file_contains_regex() {
  grep -q -E "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains_regex() {
  grep -q -E "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
#############################################

echo "> Prepare test data"

# See structure of reads:
# - https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/top_steps.html
# - https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/steps_cell_label.html
# - https://scomix.bd.com/hc/en-us/articles/360057714812-All-FAQ
# R1 is Cell Label + UMI + PolyT -> 60 bp
#   actually, CLS1 + "GTGA" + CLS2 + "GACA" + CLS3 + UMI
# R2 is the actual read -> 42 bp

# Example R1
# CLS1       Link CLS2      Link CLS3       UMI
# AAAATCCTGT GTGA AACCAAAGT GACA GATAGAGGAG CGCATGTTTATAAC


cat > reads_R1.fastq <<'EOF'
@A00226:970:H5FGVDMXY:1:1101:2645:1000 1:N:0:CAGAGAGG
AAAATCCTGTGTGAAACCAAAGTGACAGATAGAGGAGCGCATGTTTATAAC
+
FFFF:F,FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF,F,,,
@A00226:970:H5FGVDMXY:1:1101:2862:1000 1:N:0:CAGAGATG
AATAAGTGCGTGAAGAATGGAGGACAACAACTAGAATATTATGTTTGTAAA
+
:FFFF:F:FF,,FFFFF:F:F,:FFFF:FFFFF:,:,FF,FFFFFF:,FFF
EOF

# GGATGCCCAGCTAGTTTGAATTTTAGATAAACAACGAATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATTATTGGTTGTTTATCTGAGATTCAGAATTAAGCATTTTAT

# note: probably need to reverse it
cat > reads_R2.fastq <<'EOF'
@A00226:970:H5FGVDMXY:1:1101:2645:1000 2:N:0:CAGAGAGG
TCTTCGCCCGGCCAGGAATCACAAGCTCCGGGTGGATAAGGCAGCTGCTGCAGCAGCGGCACTACAAGCCA
+
FFF:FFF:F:FFFFFFFF,FFFFF:FFF:FFFFFFFFFF,FFFFFFFFFFFFFFF:FFFFFFFFFFF,FFF
@A00226:970:H5FGVDMXY:1:1101:2862:1000 2:N:0:CAGAGATG
AATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATTATTGGTT
+
F:FFFF:F,:FFFF,F:FF:F:FFFFFFFF,FF,:FFFFFFFF:FF,,F::FF::FFFFF:F:FFFFF:,F
EOF



echo "> Generate index"
CWL_FILE="$meta_resources_dir/bd_rhapsody_make_reference.cwl"
CONFIG_FILE="reference_config.yml"

cat > $CONFIG_FILE <<EOF
Genome_fasta:
  - class: File
    location: $meta_resources_dir/test_data/reference_small.fa
Gtf:
  - class: File
    location: $meta_resources_dir/test_data/reference_small.gtf
EOF

cwl-runner \
  --no-container \
  --preserve-entire-environment \
  --outdir index \
  $CWL_FILE \
  --Genome_fasta "$meta_resources_dir/test_data/reference_small.fa" \
  --Gtf "$meta_resources_dir/test_data/reference_small.gtf"
  
#########################################################################################

echo ">> Run $meta_name"
"$meta_executable" \
  --reads reads_R1.fastq \
  --reads reads_R2.fastq \
  --reference_archive "$meta_resources_dir/test_data/reference_small.tar.gz" \
  --output output \
  ${meta_cpus:+---cpus $meta_cpus} \
  ${meta_memory_mb:+---memory ${meta_memory_mb}MB}

# echo ">> Check if output exists"
# assert_file_exists "output.bam"
# assert_file_exists "log.txt"
# assert_file_exists "unmapped_r1.bam"
# assert_file_exists "unmapped_r2.bam"

# echo ">> Check if output contents are not empty"
# assert_file_not_empty "output.bam"
# assert_file_not_empty "log.txt"
# assert_file_not_empty "unmapped_r1.bam"
# assert_file_not_empty "unmapped_r2.bam"

# echo ">> Check if output contents are correct"
# assert_file_contains "log.txt" "Number of input reads \\|	2"
# assert_file_contains "log.txt" "Number of reads unmapped: too short \\|	1"
# assert_file_contains "log.txt" "Uniquely mapped reads number \\|	1"

#########################################################################################

echo "> Test successful"
