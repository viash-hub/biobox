#!/bin/bash

## VIASH START

par_input=src/multiqc/test_data/
par_output_report=woowoo.html

## VIASH END

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

# handle outputs
[[ -z "$par_output_report" ]] && no_report=true
[[ -z "$par_output_data" ]] && no_data_dir=true
[[ ! -z "$par_output_data" ]] && data_dir=true
[[ ! -z "$par_output_plots" ]] && export=true

# handle multiples
if [[ -n "$par_input" ]]; then
    IFS="," read -ra inputs <<< $par_input
    unset IFS
fi

if [[ -z "$par_include_modules" ]]; then
    include_modules=""
    IFS="," read -ra incl_modules <<< $par_include_modules
    for i in "${incl_modules[@]}"; do
        include_modules+="--include-modules $i "
    done
    unset IFS
fi

if [[ -z "$par_exclude_modules" ]]; then
    exclude_modules=""
    IFS="," read -ra excl_modules <<< $par_include_modules
    for i in "${excl_modules[@]}"; do
        exclude_modules+="--exclude-modules $i "
    done
    unset IFS
fi

if [[ -z "$par_ignore_analysis" ]]; then
    ignore=""
    IFS="," read -ra ignore_analysis <<< $par_ignore_analysis
    for i in "${ignore_analysis[@]}"; do
        ignore+="--ignore $i "
    done
    unset IFS
fi

if [[ -z "$par_ignore_samples" ]]; then
    ignore_samples=""
    IFS="," read -ra ign_samples <<< $par_ignore_samples
    for i in "${ign_samples[@]}"; do
        ignore_samples+="--ignore-samples $i "
    done
    unset IFS
fi

# run multiqc
multiqc \
    ${par_output_report:+--filename "$report_name"} \
    ${out_dir:+--outdir "$out_dir"} \
    ${no_report:+--no-report} \
    ${no_data_dir:+--no-data-dir} \
    ${data_dir:+--data-dir} \
    ${export:+--export} \
    ${par_title:+--title "$par_title"} \
    ${par_comment:+--comment "$par_comment"} \
    ${par_template:+--template "$par_template"} \
    ${par_sample_names:+--sample-names "$par_sample_names"} \
    ${par_sample_filters:+--sample-filters "$par_sample_filters"} \
    ${par_custom_css_file:+--custom-css-file "$par_custom_css_file"} \
    ${par_profile_runtime:+--profile-runtime} \
    ${par_dirs:+--dirs} \
    ${par_dirs_depth:+--dirs-depth "$par_dirs_depth"} \
    ${par_full_names:+--full-names} \
    ${par_fn_as_s_name:+--fn-as-s-name} \
    ${par_ignore_names:+--ignore-names "$par_ignore_names"} \
    ${par_ignore_symlinks:+--ignore-symlinks} \
    ${ignore_samples:+"$ignore_samples"} \
    ${ignore:+"$ignore"} \
    ${exclude_modules:+"$exclude_modules"} \
    ${include_modules:+"$include_modules"} \
    ${par_include_modules:+--include-modules "$par_include_modules"} \
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

echo "${out_dir}"
echo "${out_dir}/${report_name}_data"
echo "$par_output_data"

if [[ -n "$par_output_data" ]]; then
    mv "${out_dir}/${report_name}_data" "$par_output_data"
fi

if [[ -n "$par_output_plots" ]]; then
    mv "${out_dir}/${report_name}_plots" "$par_output_plots"
fi