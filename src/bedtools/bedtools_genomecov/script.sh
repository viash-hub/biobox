#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

[[ "$par_input_bam" == "false" ]] && unset par_input_bam
[[ "$par_depth" == "false" ]] && unset par_depth
[[ "$par_depth_zero" == "false" ]] && unset par_depth_zero
[[ "$par_bed_graph" == "false" ]] && unset par_bed_graph
[[ "$par_bed_graph_zero_coverage" == "false" ]] && unset par_bed_graph_zero_coverage
[[ "$par_split" == "false" ]] && unset par_split
[[ "$par_ignore_deletion" == "false" ]] && unset par_ignore_deletion
[[ "$par_pair_end_coverage" == "false" ]] && unset par_pair_end_coverage
[[ "$par_fragment_size" == "false" ]] && unset par_fragment_size
[[ "$par_du" == "false" ]] && unset par_du
[[ "$par_five_prime" == "false" ]] && unset par_five_prime
[[ "$par_three_prime" == "false" ]] && unset par_three_prime
[[ "$par_trackline" == "false" ]] && unset par_trackline

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
    