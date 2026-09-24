#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_bam" == "false" ]] && unset par_bam
[[ "$par_cigar_paf" == "false" ]] && unset par_cigar_paf
[[ "$par_cigar_bam" == "false" ]] && unset par_cigar_bam

# --- Argument checks -------------------------------------------------------
# --output_index is deliberately not rejected without --bam: the Nextflow
# runner fills in a default path for every output file argument, so rejecting
# it would make every PAF run fail there. Without --bam it is simply ignored.

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
  # -O bam is required: samtools sort otherwise picks the format from the
  # output extension, so --output alignment.sam would silently yield plain SAM
  # and samtools index would then fail.
  sort_args=(
    -O bam
    ${meta_cpus:+-@ "$meta_cpus"}
    -o "$par_output"
  )

  # -m is per sorting thread, and "-@ N" means N additional threads, so samtools
  # may use (N + 1) * -m in total; without -m it defaults to 768M/thread and
  # overruns the task's allocation outright. Only half the allocation is handed
  # to sort: minimap2 runs at the same time on the other side of the pipe,
  # holds the reference index in RAM and cannot spill, while sort short of
  # memory just writes more temp files.
  if [[ -n "${meta_memory_mb:-}" ]]; then
    mem_per_thread=$(( meta_memory_mb / 2 / ( ${meta_cpus:-1} + 1 ) ))
    if [[ "$mem_per_thread" -lt 128 ]]; then
      mem_per_thread=128
    fi
    sort_args+=( -m "${mem_per_thread}M" )
  fi

  echo "Running minimap2 and producing sorted BAM..."
  # -a is required for SAM/BAM output
  minimap2 \
    "${cmd_args[@]}" \
    -a \
    "$par_reference" \
    "$par_query" | \
    samtools sort "${sort_args[@]}" -

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
