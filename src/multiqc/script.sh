#!/bin/bash

# disable flags
[[ "$par_ignore_symlinks" == "false" ]] && unset ignore_symlinks
[[ "$par_dirs" == "false" ]] && unset par_dirs
[[ "$par_full_names" == "false" ]] && unset par_full_names
[[ "$par_fn_as_s_name" == "false" ]] && unset par_fn_as_s_name
[[ "$par_profile_runtime" == "false" ]] && unset par_profile_runtime
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_quiet" == "false" ]] && unset par_quiet
[[ "$par_strict" == "false" ]] && unset par_strict
[[ "$par_development" == "false" ]] && unset par_development
[[ "$par_require_logs" == "false" ]] && unset par_require_logs
[[ "$par_no_megaqc_upload" == "false" ]] && unset par_no_megaqc_upload
[[ "$par_no_ansi" == "false" ]] && unset par_no_ansi
[[ "$par_flat" == "false" ]] && unset par_flat
[[ "$par_interactive" == "false" ]] && unset par_interactive
[[ "$par_static_plot_export" == "false" ]] && unset par_static_plot_export
[[ "$par_data_dir" == "false" ]] && unset par_data_dir
[[ "$par_no_data_dir" == "false" ]] && unset par_no_data_dir
[[ "$par_zip_data_dir" == "false" ]] && unset par_zip_data_dir
[[ "$par_pdf" == "false" ]] && unset par_pdf


# handle inputs
out_dir=$(dirname "$par_output_report")
output_report_file=$(basename "$par_output_report")
report_name="${output_report_file%.*}"

# echo ============
# echo $out_dir
# echo $output_report_file
# echo $report_name
# echo ============

# handle outputs
[[ -z "$par_output_report" ]] && no_report=true
[[ -z "$par_output_data" ]] && no_data_dir=true
[[ ! -z "$par_output_data" ]] && data_dir=true
[[ ! -z "$par_output_plots" ]] && export=true

# echo ===========
# echo par_output_data: $par_output_data
# echo data_dir: $data_dir
# echo no_data_dir: $no_data_dir
# echo par_output_report: $par_output_report
# echo no_report: $no_report
# echo ===========

# allow for multiple inputs
if [[ -n "$par_input" ]]; then
    IFS=":" read -ra inputs <<< $par_input
    unset IFS
fi


# run multiqc
multiqc \
    ${par_output_report:+--filename "$report_name"} \
    ${out_dir:+--outdir "$out_dir"} \
    ${no_report:+--no-report} \
    ${export:+--export} \
    ${par_data_format:+--data-format "$par_data_format"} \
    ${par_zip_data_dir:+--zip-data-dir} \
    ${par_pdf:+--pdf} \
    ${par_interactive:+--interactive} \
    ${par_flat:+--flat} \
    ${par_verbose:+--verbose} \
    ${par_quiet:+--quiet} \
    ${par_strict:+--strict} \
    ${par_no_megaqc_upload:+--no-megaqc-upload} \
    ${par_no_ansi:+--no-ansi} \
    ${par_profile_runtime:+--profile-runtime} \
    ${par_require_logs:+--require-logs} \
    ${par_development:+--development} \
    --force \
    "${inputs[@]}"

# handle output files

if [[ -n "$par_output_data" ]]; then
    mv "${out_dir}/${report_name}_data" "$par_output_data"
fi

if [[ -n "$par_output_plots" ]]; then
    mv "${out_dir}/${report_name}_plots" "$par_output_plots"
fi