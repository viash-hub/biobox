#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

unset_if_false=(
    par_filter_control_reads
    par_filter_failing_quality
    par_skip_read_name_check
    par_sample_barcode_in_fastq_header
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Create arrays for inputs that contain multiple arguments
IFS=";" read -ra fastqs <<< "$par_fastqs"
IFS=";" read -ra read_structures <<< "$par_read_structures"
IFS=";" read -ra lane <<< "$par_lane"
IFS=";" read -ra quality_mask_threashold <<< "$par_quality_mask_threshold"
IFS=";" read -ra output_types <<< "$par_output_types"

echo "> Creating temporary directory"
# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT
echo "> Temporary directory '$TMPDIR' created"

if [ "$par_most_unmatched_to_output" -eq "0" ] && [ ! -z "$par_most_frequent_unmatched" ]; then
    echo "Requested to output 'most_frequent_unmatched' file, but 'most_unmatched_to_output' is set to 0."
    exit 1
fi

args=(
    --fastqs ${fastqs[@]}
    --sample-metadata "$par_sample_metadata"
    --output-dir "$TMPDIR"
    ${par_allowed_mismatches:+--allowed-mismatches $par_allowed_mismatches}
    ${par_min_delta:+--min-delta $par_min_delta}
    ${par_free_ns:+--free-ns $par_free_ns}
    ${par_max_no_calls:+--max-no-calls $par_max_no_calls}
    ${quality_mask_threashold:+--quality-mask-threshold "${quality_mask_threashold[*]}" }
    ${output_types:+--output-types "${output_types[*]}"}
    ${par_undetermined_sample_name:+--undetermined-sample-name ${par_undetermined_sample_name}}
    ${par_most_unmatched_to_output:+--par-most-unmatched-to-output ${par_most_unmatched_to_output}}
    ${par_override_matcher:+--override-matcher $par_override_matcher}
    ${par_metric_prefix:+--metric-prefix $par_metric_prefix}
    ${lane:+--lane "${lane[*]}"}
    ${read_structures:+--read-structures ${read_structures[*]}}
    ${par_filter_control_reads:+--filter-control-reads}
    ${par_filter_failing_quality:+--filter-failing-quality}
    ${par_skip_read_name_check:+--skip-read-name-check}
    ${par_sample_barcode_in_fastq_header:+--sample-barcode-in-fastq-header}
)

echo "> Running sgdemux with arguments: ${args[@]}"
sgdemux ${args[@]}
echo "> Done running sgdemux"

echo "> Copying FASTQ files to $par_sample_fastq"
find "$TMPDIR" -type f -name "*.fastq.gz" -exec mv '{}' "$par_sample_fastq" \;

declare -A output_files=(["metrics.tsv"]="par_metrics"
                         ["most_frequent_unmatched.tsv"]="par_most_frequent_unmatched"
                         ["sample_barcode_hop_metrics.tsv"]="par_sample_barcode_hop_metrics"
                         ["per_project_metrics.tsv"]="par_per_project_metrics"
                         ["per_sample_metrics.tsv"]="par_per_sample_metrics"
                        )

for output_file_name in "${!output_files[@]}"; do
    output_arg_variable_name=${output_files[$output_file_name]}
    destination="${!output_arg_variable_name}"
    if [ ! -z "$destination" ]; then
        echo "> Copying $output_file file to $destination"
        output_file="$TMPDIR/$output_file_name"
        if [ ! -f "$output_file" ]; then
            echo "Expected a '$output_file_name' to have been created! Exiting..."
            exit 1
        fi
        cp "$output_file" "$destination"
    fi
done

echo "> Finished!"