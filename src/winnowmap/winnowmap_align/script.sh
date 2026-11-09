#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
unset_if_false=(
  par_bam
  par_homopolymer_compressed
  par_skip_self_dual
  par_sv_off
  par_cigar_bam
  par_md_tag
  par_eqx
  par_soft_clipping
)
for par in "${unset_if_false[@]}"; do
  if [[ "${!par}" == "false" ]]; then
    unset "$par"
  fi
done

# meryl has no k to fall back on, and the k it counts with has to be the k
# winnowmap runs with: winnowmap rejects a -W list built with a different k
# ("input list of k-mers and winnowmap parameter k are inconsistent"). So when
# this component builds the list itself, the user has to say which k. No
# default is hardcoded here - that would silently drift from winnowmap's own.
if [[ -z "${par_repetitive_kmers:-}" && -z "${par_kmer_size:-}" ]]; then
  echo "Error: --kmer_size is required when --repetitive_kmers is not given:" >&2
  echo "meryl needs an explicit k to count with, and winnowmap's -k must match it." >&2
  exit 1
fi

# meta_temp_dir is a shared root (VIASH_TEMP, usually /tmp), not a per-run
# directory: fixed names under it collide between concurrent invocations, and
# nothing else clears what is written there. The meryl database runs about 3x
# the reference (measured: 130 MB for a 45 Mb reference, so several GB for a
# human genome) and is written by root under the docker engine, so leaving it
# behind fills the host with directories the caller cannot even delete. Remove
# it on every exit path, failures too.
tmp_dir=$(mktemp -d -p "$meta_temp_dir" "${meta_name}_XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

# --- Step 1: Compute repetitive k-mers with meryl (if not pre-supplied) ---
if [[ -z "${par_repetitive_kmers:-}" ]]; then
  meryl_db="$tmp_dir/merylDB"
  par_repetitive_kmers="$tmp_dir/repetitive_k${par_kmer_size}.txt"

  meryl_args=(
    count
    k="${par_kmer_size}"
    ${meta_cpus:+threads="${meta_cpus}"}
    ${meta_memory_gb:+memory="${meta_memory_gb}"}
    output "${meryl_db}"
    "${par_reference}"
  )

  echo "Computing k-mer frequencies (k=${par_kmer_size}) with meryl..."
  meryl "${meryl_args[@]}"

  # distinct=0.9998 keeps the top 0.02% most frequent k-mers. Not a value picked
  # here: it is the threshold Winnowmap's own README and paper prescribe for
  # building the -W list (https://github.com/marbl/Winnowmap#usage).
  echo "Extracting repetitive k-mers (distinct threshold = 0.9998)..."
  meryl print \
    greater-than distinct=0.9998 \
    "${meryl_db}" \
    > "${par_repetitive_kmers}"

  echo "Repetitive k-mers written to: ${par_repetitive_kmers}"
fi

# --- Step 2: Align reads with winnowmap ---
# -x goes first, as winnowmap's help asks ("preset (always applied before other
# options)"). The asm presets carry their own k (19), and -k has to keep matching
# the -W list; winnowmap re-applies an explicit -k after the preset whatever the
# order, so this is about keeping the intent readable, not a correctness fix.
#
# -k is only passed when the user asked for one. With a pre-supplied -W list and
# no --kmer_size, winnowmap uses its own default k and validates the list against
# it itself, which is one less default for this component to track.
cmd_args=(
  ${par_preset:+-x "${par_preset}"}
  -a
  -W "${par_repetitive_kmers}"
  ${par_kmer_size:+-k "${par_kmer_size}"}
  ${par_window_size:+-w "${par_window_size}"}
  ${par_homopolymer_compressed:+-H}
  ${par_split_index:+-I "${par_split_index}"}
  ${par_filter_fraction:+-f "${par_filter_fraction}"}
  ${par_max_chain_gap:+-g "${par_max_chain_gap}"}
  ${par_max_intron_length:+-G "${par_max_intron_length}"}
  ${par_max_fragment_length:+-F "${par_max_fragment_length}"}
  ${par_bandwidth:+-r "${par_bandwidth}"}
  ${par_min_minimizers:+-n "${par_min_minimizers}"}
  ${par_min_chaining_score:+-m "${par_min_chaining_score}"}
  ${par_skip_self_dual:+-X}
  ${par_secondary_ratio:+-p "${par_secondary_ratio}"}
  ${par_sv_off:+--sv-off}
  ${par_match_score:+-A "${par_match_score}"}
  ${par_mismatch_penalty:+-B "${par_mismatch_penalty}"}
  ${par_gap_open_penalty:+-O "${par_gap_open_penalty}"}
  ${par_gap_extension_penalty:+-E "${par_gap_extension_penalty}"}
  ${par_zdrop:+-z "${par_zdrop}"}
  ${par_min_peak_score:+-s "${par_min_peak_score}"}
  ${par_splice_strand:+-u "${par_splice_strand}"}
  ${par_read_group:+-R "${par_read_group}"}
  ${par_cigar_bam:+-L}
  ${par_md_tag:+--MD}
  ${par_cs_tag:+--cs="${par_cs_tag}"}
  ${par_eqx:+--eqx}
  ${par_soft_clipping:+-Y}
  ${par_minibatch_size:+-K "${par_minibatch_size}"}
  ${meta_cpus:+-t "${meta_cpus}"}
  "${par_reference}"
  "${par_query}"
)

if [[ -n "${par_bam:-}" ]]; then
  # -O bam is required: samtools sort otherwise picks the format from the
  # output extension, so --output alignment.sam would silently yield plain SAM
  # and the samtools index below would fail with "not a BGZF file".
  sort_args=(
    -O bam
    -T "$tmp_dir/sort"
    ${meta_cpus:+-@ "${meta_cpus}"}
    -o "${par_output}"
  )

  # -m is per sorting thread, and "-@ N" means N additional threads, so samtools
  # may use (N + 1) * -m in total; without -m it defaults to 768M/thread and
  # overruns the task's allocation outright. Only half the allocation is handed
  # to sort: winnowmap runs at the same time on the other side of the pipe,
  # holds the reference index in RAM and cannot spill, while sort short of
  # memory just writes more temp files.
  if [[ -n "${meta_memory_mb:-}" ]]; then
    mem_per_thread=$(( meta_memory_mb / 2 / ( ${meta_cpus:-1} + 1 ) ))
    if [[ "$mem_per_thread" -lt 128 ]]; then
      mem_per_thread=128
    fi
    sort_args+=( -m "${mem_per_thread}M" )
  fi

  echo "Running winnowmap and producing sorted BAM..."
  winnowmap "${cmd_args[@]}" | samtools sort "${sort_args[@]}" -

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
