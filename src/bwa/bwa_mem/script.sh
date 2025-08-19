#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_skip_mate_rescue" == "false" ]] && unset par_skip_mate_rescue
[[ "$par_skip_pairing" == "false" ]] && unset par_skip_pairing
[[ "$par_smart_pairing" == "false" ]] && unset par_smart_pairing
[[ "$par_ignore_alt" == "false" ]] && unset par_ignore_alt
[[ "$par_primary_5prime" == "false" ]] && unset par_primary_5prime
[[ "$par_keep_mapq" == "false" ]] && unset par_keep_mapq
[[ "$par_output_all" == "false" ]] && unset par_output_all
[[ "$par_append_comment" == "false" ]] && unset par_append_comment
[[ "$par_output_ref_header" == "false" ]] && unset par_output_ref_header
[[ "$par_soft_clipping" == "false" ]] && unset par_soft_clipping
[[ "$par_mark_secondary" == "false" ]] && unset par_mark_secondary
[[ "$par_output_xb" == "false" ]] && unset par_output_xb

# Build the command
cmd_args=(
    # Algorithm options
    ${meta_cpus:+-t "$meta_cpus"}
    ${par_min_seed_length:+-k "$par_min_seed_length"}
    ${par_band_width:+-w "$par_band_width"}
    ${par_dropoff:+-d "$par_dropoff"}
    ${par_reseed_ratio:+-r "$par_reseed_ratio"}
    ${par_seed_occurrence:+-y "$par_seed_occurrence"}
    ${par_skip_seeds:+-c "$par_skip_seeds"}
    ${par_chain_drop:+-D "$par_chain_drop"}
    ${par_seeded_bases:+-W "$par_seeded_bases"}
    ${par_mate_rescue:+-m "$par_mate_rescue"}
    ${par_skip_mate_rescue:+-S}
    ${par_skip_pairing:+-P}
    
    # Scoring options
    ${par_match_score:+-A "$par_match_score"}
    ${par_mismatch_penalty:+-B "$par_mismatch_penalty"}
    ${par_gap_open_penalty:+-O "$par_gap_open_penalty"}
    ${par_gap_extend_penalty:+-E "$par_gap_extend_penalty"}
    ${par_clipping_penalty:+-L "$par_clipping_penalty"}
    ${par_unpaired_penalty:+-U "$par_unpaired_penalty"}
    ${par_read_type:+-x "$par_read_type"}
    
    # Input/Output options
    ${par_smart_pairing:+-p}
    ${par_read_group:+-R "$par_read_group"}
    ${par_header:+-H "$par_header"}
    ${par_output:+-o "$par_output"}
    ${par_ignore_alt:+-j}
    ${par_primary_5prime:+-5}
    ${par_keep_mapq:+-q}
    ${par_batch_size:+-K "$par_batch_size"}
    ${par_verbosity:+-v "$par_verbosity"}
    ${par_min_score:+-T "$par_min_score"}
    ${par_max_hits_xa:+-h "$par_max_hits_xa"}
    ${par_score_fraction:+-z "$par_score_fraction"}
    ${par_output_all:+-a}
    ${par_append_comment:+-C}
    ${par_output_ref_header:+-V}
    ${par_soft_clipping:+-Y}
    ${par_mark_secondary:+-M}
    ${par_output_xb:+-u}
    ${par_insert_size:+-I "$par_insert_size"}
    
    # Index and input files
    "$par_index"
    "$par_reads1"
    ${par_reads2:+"$par_reads2"}
)

# Run bwa mem
bwa mem "${cmd_args[@]}"
