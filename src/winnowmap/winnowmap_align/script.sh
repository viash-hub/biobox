#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_bam" == "false" ]] && unset par_bam

# --- Argument checks -------------------------------------------------------
if [[ -z "${par_bam:-}" && -n "${par_output_index:-}" ]]; then
  echo "Error: --output_index is only valid together with --bam." >&2
  exit 1
fi

# winnowmap rejects a -W k-mer list whose k differs from its own -k
# ("input list of k-mers and winnowmap parameter k are inconsistent"), so the
# same k must drive both meryl and winnowmap.
kmer_size="${par_kmer_size:-15}"

# --- Step 1: Compute repetitive k-mers with meryl (if not pre-supplied) ---
if [[ -z "${par_repetitive_kmers:-}" ]]; then
  meryl_db="${meta_temp_dir}/merylDB"
  par_repetitive_kmers="${meta_temp_dir}/repetitive_k${kmer_size}.txt"

  meryl_args=(
    count
    k="${kmer_size}"
    ${meta_cpus:+threads="${meta_cpus}"}
    ${meta_memory_gb:+memory="${meta_memory_gb}"}
    output "${meryl_db}"
    "${par_reference}"
  )

  echo "Computing k-mer frequencies (k=${kmer_size}) with meryl..."
  meryl "${meryl_args[@]}"

  echo "Extracting repetitive k-mers (distinct threshold = 0.9998)..."
  meryl print \
    greater-than distinct=0.9998 \
    "${meryl_db}" \
    > "${par_repetitive_kmers}"

  echo "Repetitive k-mers written to: ${par_repetitive_kmers}"
fi

# --- Step 2: Align reads with winnowmap ---
cmd_args=(
  -W "${par_repetitive_kmers}"
  -k "${kmer_size}"
  -a
  ${par_preset:+-x "${par_preset}"}
  ${meta_cpus:+-t "${meta_cpus}"}
  "${par_reference}"
  "${par_query}"
)

if [[ -n "${par_bam:-}" ]]; then
  echo "Running winnowmap and producing sorted BAM..."
  winnowmap "${cmd_args[@]}" | \
    samtools sort \
      ${meta_cpus:+-@ "${meta_cpus}"} \
      -o "${par_output}" \
      -

  # Write the index to the declared output path when given, so the Nextflow
  # runner publishes it; otherwise fall back to the conventional <output>.bai.
  bam_index="${par_output_index:-${par_output}.bai}"
  samtools index \
    ${meta_cpus:+-@ "${meta_cpus}"} \
    "${par_output}" "${bam_index}"
  echo "BAM index created: ${bam_index}"
else
  echo "Running winnowmap and producing SAM..."
  winnowmap "${cmd_args[@]}" > "${par_output}"
fi

echo "Alignment finished: ${par_output}"
