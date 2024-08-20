#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset variables
unset_if_false=(
    par_input_bam
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

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

bedtools genomecov \
    ${par_depth:+-d} \
    ${par_depth_zero:+-dz} \
    ${par_bed_graph:+-bg} \
    ${par_bed_graph_zero_coverage:+-bga} \
    ${par_split:+-split} \
    ${par_ignore_deletion:+-ignoreD} \
    ${par_du:+-du} \
    ${par_five_prime:+-5} \
    ${par_three_prime:+-3} \
    ${par_trackline:+-trackline} \
    ${par_strand:+-strand "$par_strand"} \
    ${par_max:+-max "$par_max"} \
    ${par_scale:+-scale "$par_scale"} \
    ${par_trackopts:+-trackopts "$par_trackopts"} \
    ${par_input_bam:+-ibam "$par_input_bam"} \
    ${par_input:+-i "$par_input"} \
    ${par_genome:+-g "$par_genome"} \
    ${par_pair_end_coverage:+-pc} \
    ${par_fragment_size:+-fs} \
    > "$par_output"
    