#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
unset_if_false=(
    par_depth
    par_depth_zero
    par_bed_graph
    par_bed_graph_zero_coverage
    par_split
    par_ignore_deletion
    par_pair_end_coverage
    par_fragment_size
    par_du
    par_five_prime
    par_three_prime
    par_trackline
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Build command arguments array
cmd_args=()

# Input arguments (mutually exclusive)
if [[ -n "$par_input_bam" ]]; then
    cmd_args+=(-ibam "$par_input_bam")
elif [[ -n "$par_input" ]]; then
    cmd_args+=(-i "$par_input")
    # Genome file is required when using -i
    [[ -n "$par_genome" ]] && cmd_args+=(-g "$par_genome")
fi

# Output format options (mutually exclusive)
[[ -n "$par_depth" ]] && cmd_args+=(-d)
[[ -n "$par_depth_zero" ]] && cmd_args+=(-dz)
[[ -n "$par_bed_graph" ]] && cmd_args+=(-bg)
[[ -n "$par_bed_graph_zero_coverage" ]] && cmd_args+=(-bga)

# Processing options
[[ -n "$par_split" ]] && cmd_args+=(-split)
[[ -n "$par_ignore_deletion" ]] && cmd_args+=(-ignoreD)
[[ -n "$par_strand" ]] && cmd_args+=(-strand "$par_strand")
[[ -n "$par_pair_end_coverage" ]] && cmd_args+=(-pc)
[[ -n "$par_fragment_size" ]] && cmd_args+=(-fs)
[[ -n "$par_du" ]] && cmd_args+=(-du)
[[ -n "$par_five_prime" ]] && cmd_args+=(-5)
[[ -n "$par_three_prime" ]] && cmd_args+=(-3)

# Histogram options
[[ -n "$par_max" ]] && cmd_args+=(-max "$par_max")

# Scaling options
[[ -n "$par_scale" ]] && cmd_args+=(-scale "$par_scale")

# Track options
[[ -n "$par_trackline" ]] && cmd_args+=(-trackline)
if [[ -n "$par_trackopts" ]]; then
    # Handle multiple trackopts values
    IFS=";" read -ra trackopts_array <<< "$par_trackopts"
    for trackopt in "${trackopts_array[@]}"; do
        cmd_args+=(-trackopts "$trackopt")
    done
fi

# Execute bedtools genomecov
bedtools genomecov "${cmd_args[@]}" > "$par_output"
    