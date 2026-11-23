#!/bin/bash

# exit on error
set -eo pipefail

## VIASH START
## VIASH END

# Build command arguments array
cmd_args=(kraken2)

# Add database path
cmd_args+=(--db "$par_db")

# Add threads
[[ -n "$par_threads" ]] && cmd_args+=(--threads "$par_threads")

# Add confidence score
[[ -n "$par_confidence" ]] && cmd_args+=(--confidence "$par_confidence")

# Add minimum base quality
[[ -n "$par_minimum_base_quality" ]] && cmd_args+=(--minimum-base-quality "$par_minimum_base_quality")

# Add quick mode
[[ "$par_quick" == "true" ]] && cmd_args+=(--quick)

# Add minimum hit groups
[[ -n "$par_minimum_hit_groups" ]] && cmd_args+=(--minimum-hit-groups "$par_minimum_hit_groups")

# Add use-names for report
[[ "$par_use_names" == "true" ]] && cmd_args+=(--use-names)

# Add output file
cmd_args+=(--output "$par_output")

# Add report file if specified
[[ -n "$par_report" ]] && cmd_args+=(--report "$par_report")

# Add report options
[[ "$par_report_zero_counts" == "true" ]] && cmd_args+=(--report-zero-counts)
[[ "$par_report_minimizer_data" == "true" ]] && cmd_args+=(--report-minimizer-data)
[[ "$par_use_mpa_style" == "true" ]] && cmd_args+=(--use-mpa-style)

# Add classified output
[[ -n "$par_classified_out" ]] && cmd_args+=(--classified-out "$par_classified_out")

# Add unclassified output
[[ -n "$par_unclassified_out" ]] && cmd_args+=(--unclassified-out "$par_unclassified_out")

# Add paired option
[[ "$par_paired" == "true" ]] && cmd_args+=(--paired)

# Add compression flags
[[ "$par_gzip_compressed" == "true" ]] && cmd_args+=(--gzip-compressed)
[[ "$par_bzip2_compressed" == "true" ]] && cmd_args+=(--bzip2-compressed)

# Add memory mapping
[[ "$par_memory_mapping" == "true" ]] && cmd_args+=(--memory-mapping)

# Add input files
cmd_args+=("$par_reads")

# Run kraken2
"${cmd_args[@]}"
