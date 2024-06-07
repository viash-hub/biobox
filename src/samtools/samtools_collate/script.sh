#!/bin/bash

## VIASH START
## VIASH END

set -e

[[ "$par_uncompressed" == "false" ]] && unset par_uncompressed
[[ "$par_fast" == "false" ]] && unset par_fast
[[ "$par_no_pg" == "false" ]] && unset par_no_pg

samtools collate \
    "$par_input" \
    ${par_output:+-o "$par_output"} \
    ${par_reference:+-T "$par_reference"} \
    ${par_uncompressed:+-u} \
    ${par_fast:+-f} \
    ${par_working_reads:+-r "$par_working_reads"} \
    ${par_compression:+-l "$par_compression"} \
    ${par_nb_tmp_files:+-n "$par_nb_tmp_files"} \
    ${par_tmp_prefix:+-T "$par_tmp_prefix"} \
    ${par_no_pg:+-P} \
    ${par_input_fmt_option:+-O "$par_input_fmt_option"} \
    ${par_output_fmt:+-O "$par_output_fmt"} \
    ${par_output_fmt_option:+-O "$par_output_fmt_option"}
    
exit 0
