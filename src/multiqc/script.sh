#!/bin/bash

if [[ -n "$par_input" ]]; then
    IFS=":" read -ra inputs <<< $par_input
    unset IFS
fi

multiqc \
    "${inputs[@]}"

if [[ -n "$par_output_report" ]]; then
    mv multiqc_report.html "$par_output_report"
fi

# if [[ -n "$par_output_data" ]]; then
#     mv "data" "$par_output_data"
# fi


# if [[ -n "$par_output_plots" ]]; then
#     mv "plots" "$par_output_data"
# fi