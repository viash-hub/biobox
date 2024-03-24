#!/bin/bash

# Download test data from snakemake wrapper
wget https://raw.githubusercontent.com/snakemake/snakemake-wrappers/3a4f7004281efc176fd9af732ad88d00c47d432d/bio/samtools/flagstat/test/mapped/a.bam
samtools index a.bam