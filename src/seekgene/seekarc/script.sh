#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_include_introns" == "false" ]] && unset par_include_introns
[[ "$par_skip_misB" == "false" ]] && unset par_skip_misB
[[ "$par_skip_misL" == "false" ]] && unset par_skip_misL
[[ "$par_skip_multi" == "false" ]] && unset par_skip_multi
[[ "$par_skip_len" == "false" ]] && unset par_skip_len
[[ "$par_nolambda" == "false" ]] && unset par_nolambda
[[ "$par_broad" == "false" ]] && unset par_broad

# SeekArc writes results into the output directory; make sure it exists.
mkdir -p "$par_output"

# Resolve absolute paths so the tool can locate inputs regardless of its working directory.
par_rnafq1=$(realpath "$par_rnafq1")
par_rnafq2=$(realpath "$par_rnafq2")
par_atacfq1=$(realpath "$par_atacfq1")
par_atacfq2=$(realpath "$par_atacfq2")
par_refpath=$(realpath "$par_refpath")
par_output=$(realpath "$par_output")

# Build the command
cmd_args=(
  --rnafq1 "$par_rnafq1"
  --rnafq2 "$par_rnafq2"
  --atacfq1 "$par_atacfq1"
  --atacfq2 "$par_atacfq2"
  --refpath "$par_refpath"
  --samplename "$par_samplename"
  --chemistry "$par_chemistry"
  --outdir "$par_output"
  ${meta_cpus:+--core "$meta_cpus"}

  # Barcode processing options
  ${par_include_introns:+--include-introns}
  ${par_skip_misB:+--skip_misB}
  ${par_skip_misL:+--skip_misL}
  ${par_skip_multi:+--skip_multi}
  ${par_skip_len:+--skip_len}
  ${par_star_path:+--star_path "$par_star_path"}

  # Peak calling options
  ${par_qvalue:+--qvalue "$par_qvalue"}
  ${par_snapshift:+--snapshift "$par_snapshift"}
  ${par_extsize:+--extsize "$par_extsize"}
  ${par_min_len:+--min_len "$par_min_len"}
  ${par_nolambda:+--nolambda}
  ${par_broad:+--broad}
  ${par_broad_cutoff:+--broad_cutoff "$par_broad_cutoff"}

  # Cell filtering options
  ${par_min_atac_count:+--min_atac_count "$par_min_atac_count"}
  ${par_min_gex_count:+--min_gex_count "$par_min_gex_count"}
)

# Run SeekArc Tools
seekarctools arc run "${cmd_args[@]}"
