#!/bin/bash

echo ">>> Testing $meta_executable"

test_dir="${meta_resources_dir}/test_data"
gunzip -k $test_dir/index.tar.gz
tar -xf $test_dir/index.tar

echo ">>> Test 1: Paired-end reads using STAR to align reads"
"$meta_executable" \
  --id mmliver \
  --strandedness reverse \
  --paired \
  --input "${test_dir}/mmliver_1.fq.gz;${test_dir}/mmliver_2.fq.gz" \
  --index index \
  --star --star_output_genome_bam --star_gzipped_read_file --estimate_rspd --seed 1 \
  --quiet

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

# echo ">>> Test 2: Single-end reads without quality scores"
# "$meta_executable" \
#   --phred64-quals \
#   --append-names \
#   --output-genome-bam \
#   --input /$test_dir/mmliver_1.fq \
#   --index index \
#   --id mmliver_single_quals

# ####################################################################################################

# echo ">>> Test 3: Paired-end reads with quality scores, using STAR to align reads"
# "$meta_executable" \
#   --phred64-quals \
#   --fragment-length-mean 150.0 \
#   --fragment-length-sd 35.0 \
#   --output-genome-bam \
#   --calc-ci \
#   --input /$test_dir/mmliver.fq \
#   --index index \
#   --id mmliver_single_quals

echo "All tests succeeded!"
exit 0