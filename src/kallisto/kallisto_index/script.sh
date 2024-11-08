#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

unset_if_false=( par_make_unique par_aa par_distinguish )

for var in "${unset_if_false[@]}"; do
    temp_var="${!var}"
    [[ "$temp_var" == "false" ]] && unset $var
done

if [ -n "$par_kmer_size" ]; then
    if [[ "$par_kmer_size" -lt 1 || "$par_kmer_size" -gt 31 || $(( par_kmer_size % 2 )) -eq 0 ]]; then
        echo "Error: Kmer size must be an odd number between 1 and 31."
        exit 1
    fi
fi

kallisto index \
    -i "${par_index}" \
    ${par_kmer_size:+--kmer-size "${par_kmer_size}"} \
    ${par_make_unique:+--make-unique} \
    ${par_aa:+--aa} \
    ${par_distinguish:+--distinguish} \
    ${par_min_size:+--min-size "${par_min_size}"} \
    ${par_ec_max_size:+--ec-max-size "${par_ec_max_size}"} \
    ${par_d_list:+--d-list "${par_d_list}"} \
    ${meta_cpus:+--threads "${meta_cpus}"} \
    ${par_tmp:+--tmp "${par_tmp}"} \
    "${par_input}"

