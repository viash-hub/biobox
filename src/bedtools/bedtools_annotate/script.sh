#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false flags to prevent them from being passed to bedtools
unset_if_false=(
    par_counts
    par_both
    par_strand
    par_different_strand
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Build file arguments array
files_args=()
for file in "${par_files[@]}"; do
    files_args+=("$file")
done

# Build names arguments array
names_args=()
if [[ -n "${par_names}" ]]; then
    names_args+=("-names")
    for name in "${par_names[@]}"; do
        names_args+=("$name")
    done
fi

# Execute bedtools annotate
bedtools annotate \
    -i "$par_input" \
    "${names_args[@]}" \
    ${par_counts:+-counts} \
    ${par_both:+-both} \
    ${par_strand:+-s} \
    ${par_different_strand:+-S} \
    -files "${files_args[@]}" \
    > "$par_output"
