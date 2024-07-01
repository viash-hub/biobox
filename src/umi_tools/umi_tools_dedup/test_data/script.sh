#!/bin/bash

# Download test data
wget https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam
# extract 150 reads with a maximum of two reads having the same start position
samtools view -h example.bam | head -n 150 | samtools view -bS - > sample.bam
samtools index sample.bam
rm example.bam