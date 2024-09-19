#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

unset_if_false=( par_single par_single_overhang par_rf_stranded par_fr_stranded par_plaintext )

for var in "${unset_if_false[@]}"; do
    temp_var="${!var}"
    [[ "$temp_var" == "false" ]] && unset $var
done

IFS=";" read -ra input <<< $par_input

# Check if par_single is not set and ensure even number of input files
if [ -z "$par_single" ]; then
    if [ $((${#input[@]} % 2)) -ne 0 ]; then
        echo "Error: When running in paired-end mode, the number of input files must be even."
        echo "Number of input files provided: ${#input[@]}"
        exit 1
    fi
fi


mkdir -p $par_output_dir


kallisto quant \
    ${meta_cpus:+--threads $meta_cpus} \
    -i $par_index \
    ${par_gtf:+--gtf "${par_gtf}"} \
    ${par_single:+--single} \
    ${par_single_overhang:+--single-overhang} \
    ${par_fr_stranded:+--fr-stranded} \
    ${par_rf_stranded:+--rf-stranded} \
    ${par_plaintext:+--plaintext} \
    ${par_bootstrap_samples:+--bootstrap-samples "${par_bootstrap_samples}"} \
    ${par_fragment_length:+--fragment-length "${par_fragment_length}"} \
    ${par_sd:+--sd "${par_sd}"} \
    ${par_seed:+--seed "${par_seed}"} \
    -o $par_output_dir \
    ${input[*]}


