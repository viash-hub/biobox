#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_keep_chromosome" == "false" ]] && unset par_keep_chromosome
[[ "$par_chrom_first" == "false" ]] && unset par_chrom_first
[[ "$par_no_overlapping" == "false" ]] && unset par_no_overlapping
[[ "$par_bedpe_format" == "false" ]] && unset par_bedpe_format
[[ "$par_allow_beyond_chrom_end" == "false" ]] && unset par_allow_beyond_chrom_end

# Validate parameter combinations
if [[ -n "$par_exclude" ]] && [[ -n "$par_include" ]]; then
    echo "ERROR: Cannot use --exclude and --include together" >&2
    exit 1
fi

if [[ -n "$par_max_overlap" ]] && [[ -n "$par_include" ]]; then
    echo "ERROR: Cannot use --max_overlap (-f) with --include file" >&2
    exit 1
fi

# Build command arguments array
cmd_args=(
    -i "$par_input"
    -g "$par_genome"
    ${par_exclude:+-excl "$par_exclude"}
    ${par_include:+-incl "$par_include"}
    ${par_keep_chromosome:+-chrom}
    ${par_chrom_first:+-chromFirst}
    ${par_seed:+-seed "$par_seed"}
    ${par_max_overlap:+-f "$par_max_overlap"}
    ${par_no_overlapping:+-noOverlapping}
    ${par_max_tries:+-maxTries "$par_max_tries"}
    ${par_bedpe_format:+-bedpe}
    ${par_allow_beyond_chrom_end:+-allowBeyondChromEnd}
)

# Execute bedtools shuffle
bedtools shuffle "${cmd_args[@]}" > "$par_output"
