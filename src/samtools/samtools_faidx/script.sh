#!/bin/bash

## VIASH START
## VIASH END

set -e

[[ "$par_continue" == "false" ]] && unset par_continue
[[ "$par_reverse_complement" == "false" ]] && unset par_reverse_complement
[[ "$par_fastq" == "false" ]] && unset par_fastq

samtools faidx \
    "$par_fasta" \
    -o "$par_output" \
    ${par_length:+-l "$par_length"} \
    ${par_continue:+-c} \
    ${part_region_file:+-r "$par_region_file"} \
    ${par_revferse_complement:+-r} \
    ${par_mark_strand:+--mark-strand "$par_mark_strand"} \
    ${par_fai_idx:+--fai-idx "$par_fai_idx"} \
    ${par_gzi_idx:+--gzi-idx "$par_gzi_idx"} \
    ${par_fastq:+-f}

exit 0