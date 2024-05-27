#!/bin/bash

## VIASH START
## VIASH END

set -e

test_dir="${metal_executable}/test_data"


[[ "$par_error_correct_cell" == "false" ]] && unset par_error_correct_cell
[[ "$par_reconcile_pairs" == "false" ]] && unset par_reconcile_pairs
[[ "$par_3prime" == "false" ]] && unset par_3prime
[[ "$par_ignore_read_pair_suffixes" == "false" ]] && unset par_ignore_read_pair_suffixes
[[ "$par_timeit_header" == "false" ]] && unset par_timeit_header
[[ "$par_log2stderr" == "false" ]] && unset par_log2stderr

umi_tools extract \
    --stdin "$par_input" \
    --stdout "$par_output" \
    ${par_output_stats:+--output-stats "$par_output_stats"} \
    ${par_error_correct_cell:+--error-correct-cell} \
    ${par_whitelist:+--whitelist "$par_whitelist"} \
    ${par_blacklist:+--blacklist "$par_blacklist"} \
    ${par_subset_reads:+--subset-reads "$par_subset_reads"} \
    ${par_quality_filter_threshold:+--quality-filter-threshold "$par_quality_filter_threshold"} \
    ${par_quality_filter_mask:+--quality-filter-mask "$par_quality_filter_mask"} \
    ${par_quality_encoding:+--quality-encoding "$par_quality_encoding"} \
    ${par_reconcile_pairs:+--reconcile-pairs} \
    ${par_bc_pattern:+--bc-pattern "$par_bc_pattern"} \
    ${par_bc_pattern2:+--bc-pattern2 "$par_bc_pattern2"} \
    ${par_extract_method:+--extract-method "$par_extract_method"} \
    ${par_three_prime:+--3prime} \
    ${par_read2_in:+--read2-in "$par_read2_in"} \
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


exit 0