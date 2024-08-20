#!/bin/bash

wget https://github.com/pliu55/pRSEM_demo/raw/master/input/mmliver_2.fq.gz
wget https://github.com/pliu55/pRSEM_demo/raw/master/input/mmliver_1.fq.gz
wget https://github.com/pliu55/pRSEM_demo/raw/master/input/chr19.gtf.gz
wget https://github.com/pliu55/pRSEM_demo/raw/master/input/chr19.fa.gz

# decompress file without keeping the original
gunzip mmliver_1.fq.gz
gunzip mmliver_2.fq.gz

# only keep 100 reads per file
head -n 400 mmliver_1.fq > mmliver_1.fq.tmp
head -n 400 mmliver_2.fq > mmliver_2.fq.tmp

mv mmliver_1.fq.tmp mmliver_1.fq
mv mmliver_2.fq.tmp mmliver_2.fq

rm mmliver_1.fq.gz.1 mmliver_2.fq.gz.1

# prepare reference genome
./rsem-prepare-reference \
    --gtf chr19.gtf \
    chr19.fa \
    rsem_index/genome

