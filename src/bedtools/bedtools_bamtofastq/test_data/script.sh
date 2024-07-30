#!/bin/bash

# NOTE: you need to have docker running for this script to work.

# create sam file
printf "@HD\tVN:1.0\tSO:unsorted\n" > example.sam
printf "@SQ\tSN:chr1\tLN:248956422\n" >> example.sam
printf "@SQ\tSN:chr2\tLN:242193529\n" >> example.sam
printf "read1\t99\tchr1\t100\t60\t48M\t=\t150\t100\tGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAAGGTACGCGTAC\t*\tNM:i:0\tMD:Z:50\tR2:Z:TGCGTACGGTTAAGAGTACGCGTACGGTACGCGTACGCGTACGCGT\tQ2:Z:IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" >> example.sam
printf "read2\t147\tchr1\t150\t60\t48M\t=\t100\t-100\tTTTCAAAGCAGTATCGATCAAATAGTAAAGGTACGCGTACGAGTGGAA\t*\tNM:i:1\tMD:Z:10G39\tR2:Z:GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAAGGTACGCGT\tQ2:Z:IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" >> example.sam

# create bam file
# docker run -it -v `pwd`:/host biocontainers/samtools:v1.9-4-deb_cv1 bash
# cd /host
# samtools view -b example.sam > example.bam
# exit

# # create fastq files
# docker run -it -v `pwd`:/host biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bash
# cd /host

# # expected.fastq
# bedtools bamtofastq -i example.bam -fq expected.fastq

# # expected_tags.fastq
#bedtools bamtofastq -i example.bam -fq expected.fastq -tags

# # expected_fq2.fastq

# exit