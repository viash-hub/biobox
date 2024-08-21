#!/bin/bash

echo ">>> Testing $meta_executable"

test_dir="${meta_resources_dir}/test_data"

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/rsem.tar.gz
gunzip -k rsem.tar.gz
tar -xf rsem.tar

mv $test_dir/rsem $meta_resources_dir

echo ">>> Test 1: Paired-end reads using STAR to align reads"
"$meta_executable" \
    --star \
    --star_gzipped_read_file \
    --paired \
    --input "$test_dir/SRR6357070_1.fastq.gz;$test_dir/SRR6357070_2.fastq.gz" \
    --index rsem \
    --id WT_REP1 \
    --seed 1 \
    --quiet

echo ">>> Checking whether output exists"
[ ! -f "WT_REP1.genes.results" ] && echo "Gene level expression counts file does not exist!" && exit 1
[ ! -s "WT_REP1.genes.results" ] && echo "Gene level expression counts file is empty!" && exit 1
[ ! -f "WT_REP1.isoforms.results" ] && echo "Transcript level expression counts file does not exist!" && exit 1
[ ! -s "WT_REP1.isoforms.results" ] && echo "Transcript level expression counts file is empty!" && exit 1
[ ! -d "WT_REP1.stat" ] && echo "Stats file does not exist!" && exit 1

echo ">>> Check wheter output is correct"
diff $test_dir/ref.genes.results WT_REP1.genes.results || { echo "Gene level expression counts file is incorrect!"; exit 1; }
diff $test_dir/ref.isoforms.results WT_REP1.isoforms.results || { echo "Transcript level expression counts file is incorrect!"; exit 1; }
diff $test_dir/ref.cnt WT_REP1.stat/WT_REP1.cnt || { echo "Stats file is incorrect!"; exit 1; }

#####################################################################################################

echo "All tests succeeded!"
exit 0