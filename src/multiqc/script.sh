#!/bin/bash

if [[ -n "$par_input" ]]; then
    IFS="," read -ra inputs <<< $par_input
    unset IFS
    input=()
    for path in "${inputs[@]}"; do
        input+=("$path")
    done
fi


multiqc \
    "${input[@]}"
