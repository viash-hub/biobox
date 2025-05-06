#!/bin/bash

set -e

## VIASH START
## VIASH END

#########################################################################################

echo "> Prepare test data"

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

#########################################################################################

echo "> Generate index"
"$meta_executable" \
  ${meta_cpus:+---cpus $meta_cpus} \
  --index "star_index/" \
  --genome_fasta_files "genome.fasta" \
  --sjdb_gtf_file "genes.gtf" \
  --genome_sa_index_nbases 4 

files=("Genome" "Log.out" "SA" "SAindex" "chrLength.txt" "chrName.txt" "chrNameLength.txt" "chrStart.txt" "exonGeTrInfo.tab" "exonInfo.tab" "geneInfo.tab" "genomeParameters.txt" "sjdbInfo.txt" "sjdbList.fromGTF.out.tab" "sjdbList.out.tab" "transcriptInfo.tab")

echo ">> Check if output exists"
[ ! -d "star_index" ] && echo "Directory 'star_index' does not exist!" && exit 1
for file in "${files[@]}"; do
    [ ! -f "star_index/$file" ] && echo "File '$file' does not exist in 'star_index'." && exit 1
done

echo ">> Check contents of output files"
grep -q "200" "star_index/chrLength.txt" || (echo "Chromosome length in file 'chrLength.txt' is incorrect! " && exit 1)
grep -q "chr1" "star_index/chrName.txt" || (echo "Chromosome name in file 'chrName.txt' is incorrect! " && exit 1)
grep -q "chr1	200" "star_index/chrNameLength.txt" || (echo "Chromosome name in file 'chrNameLength.txt' is incorrect! " && exit 1)

echo "> Test optional parameters"
"$meta_executable" \
  ${meta_cpus:+---cpus $meta_cpus} \
  --index "star_index/" \
  --genome_fasta_files "genome.fasta" \
  --sjdb_gtf_file "genes.gtf" \
  --sjdb_overhang 100 \
  --genome_sa_index_nbases 4 \
  --sjdb_gtf_tag_exon_parent_transcript "transcript_id" \
  --sjdb_gtf_tag_exon_parent_gene "gene_id" \
  --sjdb_gtf_tag_exon_parent_gene_name "gene_name" \
  --sjdb_gtf_tag_exon_parent_gene_type "gene_type" \
  --genome_chr_bin_nbits 10 \
  --sjdb_gtf_feature_exon "exon"

files=("Genome" "Log.out" "SA" "SAindex" "chrLength.txt" "chrName.txt" "chrNameLength.txt" "chrStart.txt" "exonGeTrInfo.tab" "exonInfo.tab" "geneInfo.tab" "genomeParameters.txt" "sjdbInfo.txt" "sjdbList.fromGTF.out.tab" "sjdbList.out.tab" "transcriptInfo.tab")

echo ">> Check if output exists"
[ ! -d "star_index" ] && echo "Directory 'star_index' does not exist!" && exit 1
for file in "${files[@]}"; do
    [ ! -f "star_index/$file" ] && echo "File '$file' does not exist in 'star_index'." && exit 1
done

echo ">>> Test finished successfully"
exit 0
