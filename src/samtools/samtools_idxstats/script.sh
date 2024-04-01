#!/bin/bash

## VIASH START
## VIASH END

set -e

samtools idxstats "$par_bam" > "$par_output"