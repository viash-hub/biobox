#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

test_dir="${metal_executable}/test_data"

[[ "$par_paired" == "false" ]] && unset par_paired
[[ "$par_error_correct_cell" == "false" ]] && unset par_error_correct_cell
[[ "$par_reconcile_pairs" == "false" ]] && unset par_reconcile_pairs
[[ "$par_3prime" == "false" ]] && unset par_3prime
[[ "$par_ignore_read_pair_suffixes" == "false" ]] && unset par_ignore_read_pair_suffixes
[[ "$par_timeit_header" == "false" ]] && unset par_timeit_header
[[ "$par_log2stderr" == "false" ]] && unset par_log2stderr


function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT 

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

IFS="," read -ra input <<< "$par_input"
IFS="," read -ra pattern <<< "$par_bc_pattern"

read_count="${#input[@]}"
pattern_count="${#pattern[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - Reads: $read_count bc_patterns: $pattern_count"
    if [ "$read_count" -ne 2 ] || [ "$pattern_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        read2="$(basename -- ${input[1]})"
        umi_tools extract \
            -I "${input[0]}" \
            --read2-in "${input[1]}" \
            -S "$tmpdir/$read1" \
            --read2-out="$tmpdir/$read2" \
            ${par_umitools_extract_method:+--extract-method "$par_umitools_extract_method"} \
            --bc-pattern "${pattern[0]}" \
            --bc-pattern2 "${pattern[1]}" \
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
            cp $tmpdir/$read2 $par_fastq_2
        elif [ $par_umi_discard_read == 2 ]; then
            # discard read 2
            cp $tmpdir/$read1 $par_fastq_1
        else
            cp $tmpdir/$read1 $par_fastq_1
            cp $tmpdir/$read2 $par_fastq_2
        fi
    fi
else
    echo "Not Paired - $read_count"
    if [ "$read_count" -ne 1 ] || [ "$pattern_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        umi_tools extract \
            -I "${input[0]}" \
            -S "$tmpdir/$read1" \
            --bc-pattern "${pattern[0]}" \
            ${par_umitools_extract_method:+--extract-method "$par_umitools_extract_method"} \
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

        cp $tmpdir/$read1 $par_fastq_1
    fi
fi
