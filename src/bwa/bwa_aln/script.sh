#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_log_gap_penalty" == "false" ]] && unset par_log_gap_penalty
[[ "$par_non_iterative" == "false" ]] && unset par_non_iterative
[[ "$par_illumina_13_format" == "false" ]] && unset par_illumina_13_format
[[ "$par_input_bam" == "false" ]] && unset par_input_bam
[[ "$par_single_end_only" == "false" ]] && unset par_single_end_only
[[ "$par_use_first_read" == "false" ]] && unset par_use_first_read
[[ "$par_use_second_read" == "false" ]] && unset par_use_second_read
[[ "$par_filter_casava" == "false" ]] && unset par_filter_casava

# Resolve the index base name from the index directory. BWA takes a prefix that
# its five index files share, which is not itself a path that can be mounted or
# staged, so the directory holding them is passed in instead.
index_dir="${par_index%/}"

if [[ ! -d "$index_dir" ]]; then
    echo "Error: --index must be a directory containing the BWA index files," >&2
    echo "but '$index_dir' is not a directory." >&2
    exit 1
fi

mapfile -d '' -t bwt_files < <(find -L "$index_dir" -maxdepth 1 -name '*.bwt' -type f -print0)

if [[ "${#bwt_files[@]}" -eq 0 ]]; then
    echo "Error: No BWA index found in '$index_dir': expected a '.bwt' file." >&2
    exit 1
elif [[ "${#bwt_files[@]}" -gt 1 ]]; then
    echo "Error: Multiple BWA indices found in '$index_dir': ${bwt_files[*]}." >&2
    echo "The index directory must contain exactly one index." >&2
    exit 1
fi

index_prefix="${bwt_files[0]%.bwt}"

# Build the command
cmd_args=(
    # Algorithm options
    ${meta_cpus:+-t "$meta_cpus"}
    ${par_max_diff:+-n "$par_max_diff"}
    ${par_max_gap_opens:+-o "$par_max_gap_opens"}
    ${par_max_gap_extensions:+-e "$par_max_gap_extensions"}
    ${par_indel_end_skip:+-i "$par_indel_end_skip"}
    ${par_max_long_deletion_extensions:+-d "$par_max_long_deletion_extensions"}
    ${par_seed_length:+-l "$par_seed_length"}
    ${par_max_seed_diff:+-k "$par_max_seed_diff"}
    ${par_max_queue_entries:+-m "$par_max_queue_entries"}
    
    # Scoring options
    ${par_mismatch_penalty:+-M "$par_mismatch_penalty"}
    ${par_gap_open_penalty:+-O "$par_gap_open_penalty"}
    ${par_gap_extension_penalty:+-E "$par_gap_extension_penalty"}
    ${par_stop_search_threshold:+-R "$par_stop_search_threshold"}
    ${par_quality_threshold:+-q "$par_quality_threshold"}
    
    # Input/Output options
    ${par_output:+-f "$par_output"}
    ${par_barcode_length:+-B "$par_barcode_length"}
    ${par_log_gap_penalty:+-L}
    ${par_non_iterative:+-N}
    ${par_illumina_13_format:+-I}
    ${par_input_bam:+-b}
    ${par_single_end_only:+-0}
    ${par_use_first_read:+-1}
    ${par_use_second_read:+-2}
    ${par_filter_casava:+-Y}
    
    # Index and input file
    "$index_prefix"
    "$par_reads"
)

# Run bwa aln
bwa aln "${cmd_args[@]}"
