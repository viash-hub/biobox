#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean flags
unset_if_false=(
    par_paf_output
    par_hic
    par_meth
    par_chain_only
    par_skip_pairing
    par_no_unmapped
    par_copy_comments
    par_soft_clipping
    par_primary_5prime
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Resolve the index prefix from the index directory
index_prefix="$(basename "$(ls "${par_index%/}"/*.l2b)" .l2b)"
index_path="${par_index%/}/$index_prefix"

# Build the command
cmd_args=(
    # Common
    ${par_paf_output:+-f}
    ${meta_cpus:+-t "$meta_cpus"}
    ${par_short_read_threshold:+-l "$par_short_read_threshold"}
    ${par_read_group:+-R "$par_read_group"}
    ${par_base_tag:+-b "$par_base_tag"}
    ${par_hic:+--hic}
    ${par_meth:+--meth}

    # Mapping options
    ${par_min_seed_length:+-k "$par_min_seed_length"}
    ${par_max_seed_occurrences:+-c "$par_max_seed_occurrences"}
    ${par_max_gap_size:+-g "$par_max_gap_size"}
    ${par_bandwidth:+-w "$par_bandwidth"}
    ${par_long_bandwidth:+-W "$par_long_bandwidth"}
    ${par_min_chaining_score:+-m "$par_min_chaining_score"}
    ${par_secondary_ratio:+-p "$par_secondary_ratio"}
    ${par_max_secondary:+-N "$par_max_secondary"}
    ${par_chain_only:+--chain-only}
    ${par_mode:+-x "$par_mode"}

    # Alignment options
    ${par_match_score:+-A "$par_match_score"}
    ${par_mismatch_penalty:+-B "$par_mismatch_penalty"}
    ${par_gap_open_penalty:+-O "$par_gap_open_penalty"}
    ${par_gap_extend_penalty:+-E "$par_gap_extend_penalty"}
    ${par_min_score:+-s "$par_min_score"}

    # Paired-end options
    ${par_skip_pairing:+-P}
    ${par_rescue:+--rescue "$par_rescue"}

    # Input/Output options
    ${par_output:+-o "$par_output"}
    ${par_no_unmapped:+-u}
    ${par_secondary_output_limit:+--outn "$par_secondary_output_limit"}
    ${par_xa_max_hits:+--xa "$par_xa_max_hits"}
    ${par_copy_comments:+-y}
    ${par_soft_clipping:+-Y}
    ${par_header:+-H "$par_header"}
    ${par_primary_5prime:+-5}
    ${par_batch_size:+-K "$par_batch_size"}

    # Index and input files
    "$index_path"
    "$par_reads1"
    ${par_reads2:+"$par_reads2"}
)

# Run minibwa map
minibwa map "${cmd_args[@]}"
