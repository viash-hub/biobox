#!/bin/bash

set -e

## VIASH START
meta_executable="target/docker/star/star_align_reads/star_align_reads"
meta_resources_dir="src/star/star_align_reads"
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

cat > reads_R1.fastq <<'EOF'
@SEQ_ID1
ACGCTGCCTCATAAGCCTCACACAT
+
IIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
ACCCGCAAGATTAGGCTCCGTACAC
+
!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > reads_R2.fastq <<'EOF'
@SEQ_ID1
ATGTGTGAGGCTTATGAGGCAGCGT
+
IIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
GTGTACGGAGCCTAATCTTGCAGGG
+
!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > genome.fasta <<'EOF'
>chr1
TGGCATGAGCCAACGAACGCTGCCTCATAAGCCTCACACATCCGCGCCTATGTTGTGACTCTCTGTGAGCGTTCGTGGG
GCTCGTCACCACTATGGTTGGCCGGTTAGTAGTGTGACTCCTGGTTTTCTGGAGCTTCTTTAAACCGTAGTCCAGTCAA
TGCGAATGGCACTTCACGACGGACTGTCCTTAGCTCAGGGGA
EOF

cat > genes.gtf <<'EOF'
chr1    example_source  gene    0    50   .   +   .   gene_id "gene1"; transcript_id "transcript1";
chr1    example_source  exon    20   40   .   +   .   gene_id "gene1"; transcript_id "transcript1"; 
EOF

echo "> Generate index"
STAR \
  ${meta_cpus:+--runThreadN $meta_cpus} \
  --runMode genomeGenerate \
  --genomeDir "index/" \
  --genomeFastaFiles "genome.fasta" \
  --sjdbGTFfile "genes.gtf" \
  --genomeSAindexNbases 2
  
#########################################################################################

mkdir star_align_reads_se
cd star_align_reads_se

echo "> Run star_align_reads on SE"
"$meta_executable" \
  --input "../reads_R1.fastq" \
  --genome_dir "../index/" \
  --aligned_reads "output.sam" \
  --log "log.txt" \
  --out_reads_unmapped "Fastx" \
  --unmapped "unmapped.sam" \
  --quant_mode "TranscriptomeSAM;GeneCounts" \
  --reads_per_gene "reads_per_gene.tsv" \
  --out_sj_type Standard \
  --splice_junctions "splice_junctions.tsv" \
  --reads_aligned_to_transcriptome "transcriptome_aligned.bam" \
  ${meta_cpus:+---cpus $meta_cpus}

# TODO: Test data doesn't contain any chimeric reads yet
# --chimOutType "Junctions" \
# --chimeric_junctions "chimeric_junctions.tsv" \

echo ">> Check if output exists"
assert_file_exists "output.sam"
assert_file_exists "log.txt"
assert_file_exists "reads_per_gene.tsv"
# assert_file_exists "chimeric_junctions.tsv"
assert_file_exists "splice_junctions.tsv"
assert_file_exists "unmapped.sam"
assert_file_exists "transcriptome_aligned.bam"

echo ">> Check if output contents are not empty"
assert_file_not_empty "output.sam"
assert_file_not_empty "log.txt"
assert_file_not_empty "reads_per_gene.tsv"
# assert_file_not_empty "chimeric_junctions.tsv"
# assert_file_not_empty "splice_junctions.tsv" # TODO: test data doesn't contain any splice junctions yet
assert_file_not_empty "unmapped.sam"
assert_file_not_empty "transcriptome_aligned.bam"

echo ">> Check if output contents are correct"
assert_file_contains "log.txt" "Number of input reads \\|	2"
assert_file_contains "log.txt" "Number of reads unmapped: too short \\|	1"
assert_file_contains "log.txt" "Uniquely mapped reads number \\|	1"
assert_file_contains "reads_per_gene.tsv" "gene1	1	1	0"
assert_file_contains "reads_per_gene.tsv" "N_unmapped	1	1	1"
assert_file_contains "output.sam" "SEQ_ID1	0	chr1	17	255	25M	\\*	0	0	ACGCTGCCTCATAAGCCTCACACAT	IIIIIIIIIIIIIIIIIIIIIIIII	NH:i:1	HI:i:1	AS:i:24	nM:i:0"
assert_file_contains "unmapped.sam" "@SEQ_ID2 0:N:"
assert_file_contains "unmapped.sam" "ACCCGCAAGATTAGGCTCCGTACAC"

cd ..

#########################################################################################

mkdir star_align_reads_pe_minimal
cd star_align_reads_pe_minimal

echo ">> Run star_align_reads on PE"
"$meta_executable" \
  --input ../reads_R1.fastq \
  --input_r2 ../reads_R2.fastq \
  --genome_dir ../index/ \
  --aligned_reads output.bam \
  --log log.txt \
  --out_reads_unmapped Fastx \
  --unmapped unmapped_r1.bam \
  --unmapped_r2 unmapped_r2.bam \
  ${meta_cpus:+---cpus $meta_cpus}

echo ">> Check if output exists"
assert_file_exists "output.bam"
assert_file_exists "log.txt"
assert_file_exists "unmapped_r1.bam"
assert_file_exists "unmapped_r2.bam"

echo ">> Check if output contents are not empty"
assert_file_not_empty "output.bam"
assert_file_not_empty "log.txt"
assert_file_not_empty "unmapped_r1.bam"
assert_file_not_empty "unmapped_r2.bam"

echo ">> Check if output contents are correct"
assert_file_contains "log.txt" "Number of input reads \\|	2"
assert_file_contains "log.txt" "Number of reads unmapped: too short \\|	1"
assert_file_contains "log.txt" "Uniquely mapped reads number \\|	1"

cd ..

#########################################################################################

echo "> Test successful"
