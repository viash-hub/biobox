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

#########################################################################################

# generate index. this is not part of the unit test, but data we need to run
# the bd rhapsody pipeline.
echo "> Generate index"
CWL_FILE="$meta_resources_dir/bd_rhapsody_make_reference.cwl"
CONFIG_FILE="reference_config.yml"
REFERENCE_FILE="index/Rhap_reference.tar.gz"

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
  --outdir $(dirname "$REFERENCE_FILE") \
  "$CWL_FILE" \
  --Genome_fasta "$meta_resources_dir/test_data/reference_small.fa" \
  --Gtf "$meta_resources_dir/test_data/reference_small.gtf" \
  --Extra_STAR_params "--genomeSAindexNbases 4"


#############################################

echo "> Prepare WTA test data"

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

gzip -c > WTAreads_R1.fq.gz <<'EOF'
@A00226:970:H5FGVDMXY:1:1101:2645:1000 1:N:0:CAGAGAGG
AAAATCCTGTGTGAAACCAAAGTGACAGATAGAGGAGCGCATGTTTATAAC
+
FFFF:F,FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF,F,,,
@A00226:970:H5FGVDMXY:1:1101:2862:1000 1:N:0:CAGAGATG
AATAAGTGCGTGAAGAATGGAGGACAACAACTAGAATATTATGTTTGTAAA
+
:FFFF:F:FF,,FFFFF:F:F,:FFFF:FFFFF:,:,FF,FFFFFF:,FFF
EOF

# ENST00000607096:
# GGATGCCCAGCTAGTTTGAATTTTAGATAAACAACGAATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATTATTGGTTGTTTATCTGAGATTCAGAATTAAGCATTTTAT

# note: probably need to reverse it
# orig:
# TCTTCGCCCGGCCAGGAATCACAAGCTCCGGGTGGATAAGGCAGCTGCTGCAGCAGCGGCACTACAAGCCA
# AATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATTATTGGTT
gzip -c > WTAreads_R2.fq.gz <<'EOF'
@A00226:970:H5FGVDMXY:1:1101:2645:1000 2:N:0:CAGAGAGG
AGAAGCGGGCCGGTCCTTAGTGTTCGAGGCCCACCTATTCCGTCGACGACGTCGTCGCCGTGATGTTCGGT
+
FFF:FFF:F:FFFFFFFF,FFFFF:FFF:FFFFFFFFFF,FFFFFFFFFFFFFFF:FFFFFFFFFFF,FFF
@A00226:970:H5FGVDMXY:1:1101:2862:1000 2:N:0:CAGAGATG
TTATTAAAGCATCGTATTTATACAGGGTTCGAATCAAACCCTGTATGAATACGATTTTTTGTAATAACCAA
+
F:FFFF:F,:FFFF,F:FF:F:FFFFFFFF,FF,:FFFFFFFF:FF,,F::FF::FFFFF:F:FFFFF:,F
EOF

echo "> Prepare ABC test data"
gzip -c > ABCreads_R1.fq.gz <<'EOF'
@A01604:19:HMKLYDRXY:1:1101:2211:1000 1:N:0:CGAGGCTG
CGGTCCAGGGTGAAGGCAGCTAGACAAACAACGCGTGGACTTGTTTTAAAT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F,,:
@A01604:19:HMKLYDRXY:1:1101:2428:1000 1:N:0:CGAGGCTG
AAAGTAACCCGTGAGTACATCTAGACAGTAGAAGCACAAGTTCATTAAATA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,,,F:
EOF

gzip -c > ABCreads_R2.fq.gz <<'EOF'
@A01604:19:HMKLYDRXY:1:1101:2211:1000 2:N:0:CGAGGCTG
NTTAGTGTTCCGTTTGGAGAGTAGCTAGTTGCTGTTCGTGGTCGTTTCAAAAAAAAAAAAAAAAAAAAAAA
+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF,
@A01604:19:HMKLYDRXY:1:1101:2428:1000 2:N:0:CGAGGCTG
NGTACTGCCGGGTAGTAATGTGTTCGTAGCCGGTAATAATCTTCGTGGAAAAAAAAAAAAAAAAAAAAAAA
+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
EOF

cat > bdabseq_smallpanel.fasta <<'EOF'
>IgM|IGHM|AHS0198|pAbO Catalog_940276
TTTGGAGGGTAGCTAGTTGCAGTTCGTGGTCGTTTC
>CD19:SJ25C1|CD19|AHS0030|pAbO Catalog_940004
TAGTAATGTGTTCGTAGCCGGTAATAATCTTCGTGG
>CD278|ICOS|AHS0012|pAbO Catalog_940043
ATAGTCCGCCGTAATCGTTGTGTCGCTGAAAGGGTT
EOF

echo "> Prepare SMK test data"

gzip -c > SMKreads_R1.fq.gz <<'EOF'
@A00226:970:H5FGVDMXY:1:1101:1199:1000 1:N:0:AAGAGGCA
TCAATAGACGAGGTGAAGGTTCGCTGACAAGTCTGTACGTGTTAAACACCA
+
F:FFFFFFF:FFFFFFFFFFFF,FF,FFFF:FFF:,FFFFFFFF:F,,,F,
@A00226:970:H5FGVDMXY:1:1101:2754:1000 1:N:0:AAGAGGCA
GTTGTACCTTAGTGAGCGACCACCGACAATGGGACTCTCGGACAATTATTT
+
:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF,,:FF
EOF

gzip -c > SMKreads_R2.fq.gz <<'EOF'
@A00226:970:H5FGVDMXY:1:1101:1199:1000 2:N:0:AAGAGGCA
GTTGTCAAGATGCTACCGTTCAGAGATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGGACCGG
+
F,FFF,FFFFF:FFFF:FFF,FF:FFFFFFFFF,FFFFFF,FFFFFFFFFFFFFFF,FFF,FF:F,F,FFF
@A00226:970:H5FGVDMXY:1:1101:2754:1000 2:N:0:AAGAGGCA
GTTGTCAAGATGCTACCGTTCAGAGTGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGGA
+
F:FFFFFFFFFFFFFFFF:FFFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF
EOF


#########################################################################################

echo ">> Run $meta_name"
"$meta_executable" \
  --reads WTAreads_R1.fq.gz \
  --reads WTAreads_R2.fq.gz \
  --reads ABCreads_R1.fq.gz \
  --reads ABCreads_R2.fq.gz \
  --reads SMKreads_R1.fq.gz \
  --reads SMKreads_R2.fq.gz \
  --reference_archive "$REFERENCE_FILE" \
  --abseq_reference bdabseq_smallpanel.fasta \
  --output output \
  ${meta_cpus:+---cpus $meta_cpus} \
  ${meta_memory_mb:+---memory ${meta_memory_mb}MB} \
  --cell_calling_data mRNA \
  --exact_cell_count 2 \
  --expected_cell_count 2 \
  --exclude_intronic_reads false \
  --tag_names 1-Jurkat \
  --tag_names 2-Ramos \
  --tag_names 3-THP1

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

# TODO: add test with ABC, VDJ, SMK, and ATAC

echo "> Test successful"
