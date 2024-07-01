#!/bin/bash

## VIASH START
## VIASH END

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/transcriptome.fasta

head -n 23 transcriptome.fasta > test.fasta # kepp only 4 first entries of the file for testing.

rm transcriptome.fasta