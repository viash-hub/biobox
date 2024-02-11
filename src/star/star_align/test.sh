#!/bin/bash

set -e

## VIASH START
meta_executable="target/docker/star/star_align/star_align"
meta_resources_dir="src/star/star_align"
## VIASH END


echo "> Generate index"
STAR \
  ${meta_cpus:+--runThreadN $meta_cpus} \
  --runMode genomeGenerate \
  --genomeDir "index/" \
  --genomeFastaFiles "$meta_resources_dir/test_data/genome.fasta" #\
  # --sjdbGTFfile "$meta_resources_dir/test_data/genes.gtf"
  
#########################################################################################
mkdir star_align_se
cd star_align_se

echo "> Run star_align on SE"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/a_R1.1.fastq" \
  --input "$meta_resources_dir/test_data/a_R1.2.fastq" \
  --genomeDir "../index/" \
  --aligned_reads "output.bam" \
  --log "log.txt" \
  --outReadsUnmapped "Fastx" \
  --unmapped "unmapped.bam" \
  ${meta_cpus:+---cpus $meta_cpus}

  # optional outputs
  # \
  # --quantMode "TranscriptomeSAM;GeneCounts" \
  # --reads_per_gene "reads_per_gene.tsv" \
  # --chimOutType "Junctions" \
  # --chimeric_junctions "chimeric_junctions.tsv" \
  # --outSJtype Standard \
  # --splice_junctions "splice_junctions.tsv"

# ↑ automate depending on desired outputs?

echo ">> Check if output exists"
[ ! -f "output.bam" ] && echo ">> output.bam does not exist" && exit 1
[ ! -f "log.txt" ] && echo ">> log.txt does not exist" && exit 1

# [ ! -f "reads_per_gene.tsv" ] && echo ">> reads_per_gene.tsv does not exist" && exit 1
# [ ! -f "chimeric_junctions.tsv" ] && echo ">> chimeric_junctions.tsv does not exist" && exit 1
# [ ! -f "splice_junctions.tsv" ] && echo ">> splice_junctions.tsv does not exist" && exit 1
[ ! -f "unmapped.bam" ] && echo ">> unmapped.bam does not exist" && exit 1

cd ..

#########################################################################################
mkdir star_align_pe_minimal
cd star_align_pe_minimal

echo ">> Run star_align on PE"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/a_R1.1.fastq" \
  --input "$meta_resources_dir/test_data/a_R1.2.fastq" \
  --input_r2 "$meta_resources_dir/test_data/a_R2.1.fastq" \
  --input_r2 "$meta_resources_dir/test_data/a_R2.2.fastq" \
  --genomeDir "../index/" \
  --aligned_reads "output.bam" \
  --log "log.txt" \
  --outReadsUnmapped "Fastx" \
  --unmapped "unmapped.bam" \
  --unmapped_r2 "unmapped_r2.bam" \
  ${meta_cpus:+---cpus $meta_cpus}

echo ">> Check if output exists"
[ ! -f "output.bam" ] && echo ">> output.bam does not exist" && exit 1
[ ! -f "log.txt" ] && echo ">> log.txt does not exist" && exit 1
[ ! -f "unmapped.bam" ] && echo ">> unmapped.bam does not exist" && exit 1
[ ! -f "unmapped_r2.bam" ] && echo ">> unmapped_r2.bam does not exist" && exit 1

cd ..
#########################################################################################

echo "> Test successful"