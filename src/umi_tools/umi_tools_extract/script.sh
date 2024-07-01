#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

test_dir="${metal_executable}/test_data"

[[ "$par_paired" == "false" ]] && unset par_paired
[[ "$par_error_correct_cell" == "false" ]] && unset par_error_correct_cell
[[ "$par_reconcile_pairs" == "false" ]] && unset par_reconcile_pairs
[[ "$par_three_prime" == "false" ]] && unset par_three_prime
[[ "$par_ignore_read_pair_suffixes" == "false" ]] && unset par_ignore_read_pair_suffixes
[[ "$par_timeit_header" == "false" ]] && unset par_timeit_header
[[ "$par_log2stderr" == "false" ]] && unset par_log2stderr


if [ -n "$par_paired" ]; then
    # For paired-end rendscheck that we have two read files, two patterns
    if [ -z "$par_input" ] || [ -z "$par_read2_in" ] ||
       [ -z "$par_bc_pattern" ] || [ -z "$par_bc_pattern2" ]; 
        then
        echo "Paired end input requires two read files, two UMI patterns, and two output files"
        exit 1
    fi
else 
    # For single-end reads, check that we have only one read file, one pattern
    if [ -n "$par_read2_in" ] || [ -n "$par_bc_pattern2" ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    elif [ "$par_umi_discard_read" != 0 ]; then
        echo "umi_discard_read is only valid when processing paired end reads."
        exit 1
    fi
fi

umi_tools extract \
    -I "$par_input" \
    ${par_read2_in:+ --read2-in "$par_read2_in"} \
    -S "$par_read1_out" \
    ${par_read2_out:+--read2-out "$par_read2_out"} \
    ${par_umitools_extract_method:+--extract-method "$par_umitools_extract_method"} \
    --bc-pattern "$par_bc_pattern" \
    ${par_bc_pattern2:+ --bc-pattern2 "$par_bc_pattern2"} \
    ${par_umitools_umi_separator:+--umi-separator "$par_umitools_umi_separator"} \
    ${par_output_stats:+--output-stats "$par_output_stats"} \
    ${par_error_correct_cell:+--error-correct-cell} \
    ${par_whitelist:+--whitelist "$par_whitelist"} \
    ${par_blacklist:+--blacklist "$par_blacklist"} \
    ${par_subset_reads:+--subset-reads "$par_subset_reads"} \
    ${par_quality_filter_threshold:+--quality-filter-threshold "$par_quality_filter_threshold"} \
    ${par_quality_filter_mask:+--quality-filter-mask "$par_quality_filter_mask"} \
    ${par_quality_encoding:+--quality-encoding "$par_quality_encoding"} \
    ${par_reconcile_pairs:+--reconcile-pairs} \
    ${par_three_prime:+--3prime} \
    ${par_filtered_out:+--filtered-out "$par_filtered_out"} \
    ${par_filtered_out2:+--filtered-out2 "$par_filtered_out2"} \
    ${par_ignore_read_pair_suffixes:+--ignore-read-pair-suffixes} \
    ${par_random_seed:+--random-seed "$par_random_seed"} \
    ${par_temp_dir:+--temp-dir "$par_temp_dir"} \
    ${par_compresslevel:+--compresslevel "$par_compresslevel"} \
    ${par_timeit:+--timeit "$par_timeit"} \
    ${par_timeit_name:+--timeit-name "$par_timeit_name"} \
    ${par_timeit_header:+--timeit-header} \
    ${par_log:+--log "$par_log"} \
    ${par_log2stderr:+--log2stderr} \
    ${par_verbose:+--verbose "$par_verbose"} \
    ${par_error:+--error "$par_error"}


if [ $par_umi_discard_read == 1 ]; then
    # discard read 1
    rm "$par_read1_out"
elif [ $par_umi_discard_read == 2 ]; then
    # discard read 2 (-f to bypass file existence check)
    rm -f "$par_read2_out"
fi