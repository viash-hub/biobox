#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_bam" == "false" ]] && unset par_bam
[[ "$par_cigar_paf" == "false" ]] && unset par_cigar_paf
[[ "$par_cigar_bam" == "false" ]] && unset par_cigar_bam

# --- Argument checks -------------------------------------------------------
if [[ -z "${par_bam:-}" && -n "${par_output_index:-}" ]]; then
  echo "Error: --output_index is only valid together with --bam." >&2
  exit 1
fi

# -L only affects SAM/BAM records and -c only affects PAF, so rejecting the
# wrong combination is clearer than letting minimap2 silently ignore the flag.
if [[ -z "${par_bam:-}" && -n "${par_cigar_bam:-}" ]]; then
  echo "Error: --cigar_bam is only valid together with --bam." >&2
  exit 1
fi

if [[ -n "${par_bam:-}" && -n "${par_cigar_paf:-}" ]]; then
  echo "Error: --cigar_paf is only valid without --bam; BAM records always carry a CIGAR." >&2
  exit 1
fi

# --- Align -----------------------------------------------------------------
# minimap2 recommends giving -x first, so that the preset does not override the
# options that follow it.
cmd_args=(
  ${par_preset:+-x "$par_preset"}
  ${meta_cpus:+-t "$meta_cpus"}
  ${par_cigar_paf:+-c}
  ${par_cigar_bam:+-L}
)

if [[ -n "${par_bam:-}" ]]; then
  echo "Running minimap2 and producing sorted BAM..."
  # -a is required for SAM/BAM output
  minimap2 \
    "${cmd_args[@]}" \
    -a \
    "$par_reference" \
    "$par_query" | \
    samtools sort \
      ${meta_cpus:+-@ "$meta_cpus"} \
      -o "$par_output" \
      -

  # Write the index to the declared output path when given, so the Nextflow
  # runner publishes it; otherwise fall back to the conventional <output>.bai.
  bam_index="${par_output_index:-${par_output}.bai}"
  samtools index \
    ${meta_cpus:+-@ "$meta_cpus"} \
    "$par_output" "$bam_index"
  echo "BAM index created: $bam_index"
else
  echo "Running minimap2 and producing PAF..."
  minimap2 \
    "${cmd_args[@]}" \
    "$par_reference" \
    "$par_query" \
    > "$par_output"
fi

echo "Alignment finished: $par_output"
