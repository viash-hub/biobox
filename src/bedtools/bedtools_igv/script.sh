#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_collapse_reads" == "false" ]] && unset par_collapse_reads
[[ "$par_use_name" == "false" ]] && unset par_use_name

# Build command arguments array
cmd_args=(
  -i "$par_input"
  ${par_output_path:+-path "$par_output_path"}
  ${par_session_file:+-sess "$par_session_file"}
  ${par_sort_reads:+-sort "$par_sort_reads"}
  ${par_collapse_reads:+-clps}
  ${par_use_name:+-name}
  ${par_flank_size:+-slop "$par_flank_size"}
  ${par_image_format:+-img "$par_image_format"}
)

# Execute bedtools igv
bedtools igv "${cmd_args[@]}" > "$par_output"
