#!/bin/bash

# Download test data
wget https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam
samtools view -b -o sample.bam -s 0.00005 example.bam
samtools index sample.bam > sample.bam.bai
rm example.bam