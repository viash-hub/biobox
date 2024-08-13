#!/bin/bash

## VIASH START
## VIASH END

set -exo pipefail

unset_if_false=(
   par_error_correct_cell
   par_reconcile_pairs
   par_three_prime
   par_ignore_read_pair_suffixes
   par_timeit_header
   par_log2stderr
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done


# Check if we have the correct number of input files and patterns for paired-end or single-end reads

# For paired-end rends, check that we have two read files, two patterns
# Check for paired-end inputs
if [ -n "$par_input" ] && [ -n "$par_read2_in" ]; then
    # Paired-end checks: Ensure both UMI patterns are provided
    if [ -z "$par_bc_pattern" ] || [ -z "$par_bc_pattern2" ]; then
        echo "Paired end input requires two UMI patterns."
        exit 1
    fi
elif [ -n "$par_input" ]; then
    # Single-end checks: Ensure no second read or UMI pattern for the second read is provided
    if [ -n "$par_bc_pattern2" ]; then
        echo "Single end input requires only one read file and one UMI pattern."
        exit 1
    fi
    # Check that discard_read is not set or set to 0 for single-end reads
    if [ -n "$par_umi_discard_read" ] && [ "$par_umi_discard_read" != 0 ]; then
        echo "umi_discard_read is only valid when processing paired end reads."
        exit 1
    fi
else
    # No inputs provided
    echo "No input files provided."
    exit 1
fi




umi_tools extract \
    -I "$par_input" \
    ${par_read2_in:+ --read2-in "$par_read2_in"} \
    -S "$par_output" \
    ${par_read2_out:+--read2-out "$par_read2_out"} \
    ${par_extract_method:+--extract-method "$par_extract_method"} \
    --bc-pattern "$par_bc_pattern" \
    ${par_bc_pattern2:+ --bc-pattern2 "$par_bc_pattern2"} \
    ${par_umi_separator:+--umi-separator "$par_umi_separator"} \
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


if [ "$par_umi_discard_read" == 1 ]; then
    # discard read 1
    rm "$par_read1_out"
elif [ "$par_umi_discard_read" == 2 ]; then
    # discard read 2 (-f to bypass file existence check)
    rm -f "$par_read2_out"
fi