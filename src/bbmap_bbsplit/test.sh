#!/bin/bash

echo ">>> Test $meta_functionality_name"

cat > bbsplit_fasta_list.txt << HERE
sarscov2,${meta_resources_dir}/test_data/sarscov2.fa
human,${meta_resources_dir}/test_data/human.fa
HERE

echo ">>> Building BBSplit index"
"${meta_executable}" \
  --primary_ref "${meta_resources_dir}/test_data/genome.fasta" \
  --ref_fasta_list bbsplit_fasta_list.txt \
  --only_build_index \
  --index "BBSplit_index" 

echo ">>> Check whether output exists"
[ ! -d "BBSplit_index" ] && echo "BBSplit index does not exist!" && exit 1
[ -z "$(ls -A 'BBSplit_index')" ] && echo "BBSplit index is empty!" && exit 1

echo ">>> Filtering ribosomal RNA reads"

echo ">>> Testing with single-end reads and primary/non-primary FASTA files"
"${meta_executable}" \
  --input "${meta_resources_dir}/test_data/SRR6357070_1.fastq.gz" \
  --primary_ref "${meta_resources_dir}/test_data/genome.fasta" \
  --ref_fasta_list bbsplit_fasta_list.txt \
  --fastq_1 filtered_SRR6357070_1.fastq.gz

echo ">>> Check whether output exists"
[ ! -f "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered reads file does not exist!" && exit 1
[ ! -s "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered reads file is empty!" && exit 1

rm filtered_SRR6357070_1.fastq.gz

echo ">>> Testing with paired-end reads and primary/non-primary FASTA files"
"${meta_executable}" \
  --paired \
  --input "${meta_resources_dir}/test_data/SRR6357070_1.fastq.gz,${meta_resources_dir}/test_data/SRR6357070_2.fastq.gz" \
  --primary_ref "${meta_resources_dir}/test_data/genome.fasta" \
  --ref_fasta_list "bbsplit_fasta_list.txt" \
  --fastq_1 "filtered_SRR6357070_1.fastq.gz" \
  --fastq_2 "filtered_SRR6357070_2.fastq.gz"

echo ">>> Check whether output exists"
[ ! -f "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered read 1 file does not exist!" && exit 1
[ ! -s "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered read 1 file is empty!" && exit 1
[ ! -f "filtered_SRR6357070_2.fastq.gz" ] && echo "Filtered read 2 file does not exist!" && exit 1
[ ! -s "filtered_SRR6357070_2.fastq.gz" ] && echo "Filtered read 2 file is empty!" && exit 1

rm filtered_SRR6357070_1.fastq.gz filtered_SRR6357070_2.fastq.gz

echo ">>> Testing with single-end reads and BBSplit index"
"${meta_executable}" \
  --input "${meta_resources_dir}/test_data/SRR6357070_1.fastq.gz" \
  --index "BBSplit_index" \
  --fastq_1 "filtered_SRR6357070_1.fastq.gz"

echo ">>> Check whether output exists"
[ ! -f "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered reads file does not exist!" && exit 1
[ ! -s "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered reads file is empty!" && exit 1

echo ">>> Testing with paired-end reads and BBSplit index"
"${meta_executable}" \
  --paired \
  --input "${meta_resources_dir}/test_data/SRR6357070_1.fastq.gz,${meta_resources_dir}/test_data/SRR6357070_2.fastq.gz" \
  --index "BBSplit_index" \
  --fastq_1 "filtered_SRR6357070_1.fastq.gz" \
  --fastq_2 "filtered_SRR6357070_2.fastq.gz"

echo ">>> Check whether output exists"
[ ! -f "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered read 1 file does not exist!" && exit 1
[ ! -s "filtered_SRR6357070_1.fastq.gz" ] && echo "Filtered read 1 file is empty!" && exit 1
[ ! -f "filtered_SRR6357070_2.fastq.gz" ] && echo "Filtered read 2 file does not exist!" && exit 1
[ ! -s "filtered_SRR6357070_2.fastq.gz" ] && echo "Filtered read 2 file is empty!" && exit 1

rm filtered_SRR6357070_1.fastq.gz filtered_SRR6357070_2.fastq.gz

echo "All tests succeeded!"
exit 0