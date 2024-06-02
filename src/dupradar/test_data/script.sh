#!/bin/bash

wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/genes.gtf
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam 

# subsample the bam file
samtools view -s 0.02 -b wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam > sample.bam

rm wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam