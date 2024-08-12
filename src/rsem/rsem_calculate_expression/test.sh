#!/bin/bash

echo ">>> Testing $meta_executable"

test_dir="${meta_resources_dir}/test_data"

gunzip "${test_dir}/rsem.tar.gz"
tar -xf "${test_dir}/rsem.tar"

echo ">>> Calculating expression"
"$meta_executable" \
  --id WT_REP1 \
  --strandedness reverse \
  --paired true \
  --input "${test_dir}/SRR6357070_1.fastq.gz;${test_dir}/SRR6357070_2.fastq.gz" \
  --index rsem \
  --extra_args "--star --star-output-genome-bam --star-gzipped-read-file --estimate-rspd --seed 1" \
  --logs WT_REP1.log \

echo ">>> Checking whether output exists"
[ ! -f "WT_REP1.genes.results" ] && echo "Gene level expression counts file does not exist!" && exit 1
[ ! -s "WT_REP1.genes.results" ] && echo "Gene level expression counts file is empty!" && exit 1
[ ! -f "WT_REP1.isoforms.results" ] && echo "Transcript level expression counts file does not exist!" && exit 1
[ ! -s "WT_REP1.isoforms.results" ] && echo "Transcript level expression counts file is empty!" && exit 1
[ ! -f "WT_REP1.stat" ] && echo "Stats file does not exist!" && exit 1
[ ! -s "WT_REP1.stat" ] && echo "Stats file is empty!" && exit 1
[ ! -f "WT_REP1.log" ] && echo "Log file does not exist!" && exit 1
[ ! -s "WT_REP1.log" ] && echo "Log file is empty!" && exit 1


####################################################################################################

echo ">>> Test 2: Single-end reads without quality scores"
"$meta_executable" \
  --no-qualities \
  --input "${test_dir}/SRR6357070_1.fastq.gz" \
  --index rsem \
  --id WT_REP1_no_qual

####################################################################################################

echo ">>> Test 3: Paired-end reads with quality scores, using STAR to align reads"
"$meta_executable" \
  --star \
  --star-gzipped-read-file \
  --paired-end \
  --input "${test_dir}/SRR6357070_1.fastq.gz;${test_dir}/SRR6357070_2.fastq.gz" \
  --index rsem \
  WT_REP1_star

echo "All tests succeeded!"
exit 0