#!/bin/bash

# create sam file
printf "@SQ\tSN:chr2:172936693-172938111\tLN:1418\n" > example.sam
printf "my_read\t99\tchr2:172936693-172938111\t129\t60\t100M\t=\t429\t400\tCTAACTAGCCTGGGAAAAAAGGATAGTGTCTCTCTGTTCTTTCATAGGAAATGTTGAATCAGACCCCTACTGGGAAAAGAAATTTAATGCATATCTCACT\t*\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:37\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:100\n" >> example.sam
printf "my_read\t147\tchr2:172936693-172938111\t429\t60\t100M\t=\t129\t-400\tTCGAGCTCTGCATTCATGGCTGTGTCTAAAGGGCATGTCAGCCTTTGATTCTCTCTGAGAGGTAATTATCCTTTTCCTGTCACGGAACAACAAATGATAG\t*\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:37\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:100\n" >> example.sam

# create bam file
# samtools view -b example.sam > example.bam

# create fastq files
# bedtools bamtofastq -i example.bam -fq expected.fastq
# bedtools bamtofastq -i example.bam -fq expected_1.fastq -fq2 expected_2.fastq
