#!/bin/bash

set -exo pipefail 


inner_distance.py \
    -i $par_input_file \
    -r $par_refgene \
    -o $par_output_prefix \
    ${par_sample_size:+-k "${par_sample_size}"} \
    ${par_lower_bound:+-l "${par_lower_bound}"} \
    ${par_upper_bound:+-u "${par_upper_bound}"} \
    ${par_step:+-s "${par_step}"} \
    ${par_mapq:+-q "${par_mapq}"} \
> stdout.txt

if [[ -n $par_output_stats ]]; then head -n 2 stdout.txt > $par_output_stats; fi


[[ -n "$par_output_dist" && -f "$par_output_prefix.inner_distance.txt" ]] && mv $par_output_prefix.inner_distance.txt $par_output_dist
[[ -n "$par_output_plot" && -f "$par_output_prefix.inner_distance_plot.pdf" ]] && mv $par_output_prefix.inner_distance_plot.pdf $par_output_plot
[[ -n "$par_output_plot_r" && -f "$par_output_prefix.inner_distance_plot.r" ]] && mv $par_output_prefix.inner_distance_plot.r $par_output_plot_r
[[ -n "$par_output_freq" && -f "$par_output_prefix.inner_distance_freq.txt" ]] && mv $par_output_prefix.inner_distance_freq.txt $par_output_freq

exit 0