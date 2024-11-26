#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
unset_if_false=(
    par_demux_only
    par_detect_adapters
    par_error_on_missing
    par_group_fastq
    par_legacy_fastq
    par_no_error_on_invalid
    par_no_projects
    par_qc_only
    par_split_lanes
    par_skip_qc_report
    par_strict
    par_force_index_orientation
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# NOTE: --preparation-workflow is bugged in bases2fastq
args=(
    ${par_demux_only:+--demux-only}
    ${par_detect_adapters:+--detect-adapters}
    ${par_error_on_missing:+--error-on-missing}
    ${par_group_fastq:+--group-fastq}
    ${par_legacy_fastq:+--legacy-fastq}
    ${par_no_error_on_invalid:+--no-error-on-invalid}
    ${par_no_projects:+--no-projects}
    ${par_split_lanes:+--split-lanes}
    ${par_strict:+--strict}
    ${par_force_index_orientation:+--force-index-orientation}
    ${par_chemistry_version:+--chemistry-version "$par_chemistry_version"}
    ${par_filter_mask:+--filter-mask "$par_filter_mask"}
    ${par_flowcell_id:+--flowcell-id "$par_flowcell_id"}
    ${par_i1_cycles:+--i1-cycles "$par_i1_cycles"}
    ${par_i2_cycles:+--i2-cycles "$par_i2_cycles"}
    ${par_r1_cycles:+--r1-cycles "$par_r1_cycles"}
    ${par_r2_cycles:+--r2-cycles "$par_r2_cycles"}
    ${par_kit_configuration:+--kit-configuration "$par_kit_configuration"}
    ${par_log_level:+--log-level "$par_log_level"}
    ${par_num_unassigned:+--num-unassigned "$par_num_unassigned"}
    ${par_preparation_workflow:+--preparation-workflow "$par_preparation_workflow"}
    ${meta_cpus:+--num-threads "$meta_cpus"}
    ${par_run_manifest:+--run-manifest "$par_run_manifest"}
)

# Create arrays for inputs that contain multiple arguments
IFS=";" read -ra exclude_tile <<< "$par_exclude_tile"
IFS=";" read -ra include_tile <<< "$par_include_tile"

if [ -z "$par_report" ]; then
    args+=( --skip-qc-report )
fi

for arg_value in "${exclude_tile[@]}"; do
    args+=( "--exclude-tile" "$arg_value" )
done

for arg_value in "${include_tile[@]}"; do
    args+=( "--include-tile" "$arg_value" )
done

echo "> Creating temporary directory."
# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
echo "> Created $TMPDIR"
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT

args+=( "$par_analysis_directory" "$TMPDIR")
echo "> Running bases2fastq with arguments: ${args[@]}"
bases2fastq ${args[@]}
echo "> Done running sgdemux"

echo "> Output folder:"
tree "$TMPDIR"

echo "> Moving FASTQ files into final output directory"
mkdir -p "$par_output_directory/"
mv "$TMPDIR"/Samples/* --target-directory="$par_output_directory"

if [ ! -z "$par_report" ]; then
  echo "> Moving HTML report to the output ($par_report)"
  mv "$TMPDIR"/*.html "$par_report"
else
  echo " > Leaving reports alone"
fi

# Logs is everything else
if [ ! -z "$par_logs" ]; then
  mkdir -p "$par_logs"
  echo "> Moving logs to their own location ($par_logs)"
  mv "$TMPDIR/"* "$par_logs/"
else
  echo "> Not moving logs"
fi



