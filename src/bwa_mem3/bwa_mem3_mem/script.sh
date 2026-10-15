#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_smem_dedup" == "false" ]] && unset par_smem_dedup
[[ "$par_huge_pages" == "false" ]] && unset par_huge_pages
[[ "$par_keep_contained_ext" == "false" ]] && unset par_keep_contained_ext
[[ "$par_adaptive_band" == "false" ]] && unset par_adaptive_band
[[ "$par_no_adaptive_band" == "false" ]] && unset par_no_adaptive_band
[[ "$par_no_band_cert" == "false" ]] && unset par_no_band_cert
[[ "$par_extend_csub" == "false" ]] && unset par_extend_csub
[[ "$par_rescue_skip" == "false" ]] && unset par_rescue_skip
[[ "$par_skip_mate_rescue" == "false" ]] && unset par_skip_mate_rescue
[[ "$par_skip_pairing" == "false" ]] && unset par_skip_pairing
[[ "$par_hic" == "false" ]] && unset par_hic
[[ "$par_fast" == "false" ]] && unset par_fast
[[ "$par_compat_allow_divergent" == "false" ]] && unset par_compat_allow_divergent
[[ "$par_smart_pairing" == "false" ]] && unset par_smart_pairing
[[ "$par_ignore_alt" == "false" ]] && unset par_ignore_alt
[[ "$par_primary_5prime" == "false" ]] && unset par_primary_5prime
[[ "$par_keep_mapq" == "false" ]] && unset par_keep_mapq
[[ "$par_output_xb" == "false" ]] && unset par_output_xb
[[ "$par_output_all" == "false" ]] && unset par_output_all
[[ "$par_append_comment" == "false" ]] && unset par_append_comment
[[ "$par_output_ref_header" == "false" ]] && unset par_output_ref_header
[[ "$par_soft_clipping" == "false" ]] && unset par_soft_clipping
[[ "$par_mark_secondary" == "false" ]] && unset par_mark_secondary
[[ "$par_chimera_qc" == "false" ]] && unset par_chimera_qc
[[ "$par_proper_pair_from_emitted" == "false" ]] && unset par_proper_pair_from_emitted
[[ "$par_legacy_reader" == "false" ]] && unset par_legacy_reader

# --extend-mate-concordant takes an optional argument: the bare flag means auto,
# whereas '=auto' would be parsed by atoi() as 0 (off).
if [[ "$par_extend_mate_concordant" == "auto" ]]; then
  extend_mate_concordant_arg="--extend-mate-concordant"
elif [[ -n "$par_extend_mate_concordant" ]]; then
  extend_mate_concordant_arg="--extend-mate-concordant=$par_extend_mate_concordant"
fi

# Build the command
# Options with an optional argument (--bam, --rescue-kmer, --extend-mate-concordant,
# --meth, --meth-seed-prune) only bind their value when written as --opt=VALUE.
cmd_args=(
  ${meta_cpus:+-t "$meta_cpus"}

  # Output
  ${par_output:+-o "$par_output"}

  # BAM options
  ${par_bam_compression:+--bam="$par_bam_compression"}
  ${par_bam_threads:+--bam-threads "$par_bam_threads"}

  # Algorithm options
  ${par_min_seed_length:+-k "$par_min_seed_length"}
  ${par_band_width:+-w "$par_band_width"}
  ${par_dropoff:+-d "$par_dropoff"}
  ${par_reseed_ratio:+-r "$par_reseed_ratio"}
  ${par_split_width:+-s "$par_split_width"}
  ${par_seed_occurrence:+-y "$par_seed_occurrence"}
  ${par_skip_seeds:+-c "$par_skip_seeds"}
  ${par_smem_dedup:+--smem-dedup}
  ${par_dedup:+--dedup "$par_dedup"}
  ${par_dedup_reads:+--dedup-reads "$par_dedup_reads"}
  ${par_ks_dedup:+--ks-dedup "$par_ks_dedup"}
  ${par_huge_pages:+--huge-pages}
  ${par_keep_contained_ext:+--keep-contained-ext}
  ${par_max_extend_chains:+--max-extend-chains "$par_max_extend_chains"}
  ${par_adaptive_band:+--adaptive-band}
  ${par_no_adaptive_band:+--no-adaptive-band}
  ${par_no_band_cert:+--no-band-cert}
  ${extend_mate_concordant_arg:+"$extend_mate_concordant_arg"}
  ${par_extend_tie_frac:+--extend-tie-frac "$par_extend_tie_frac"}
  ${par_extend_csub:+--extend-csub}
  ${par_extend_tie_floor:+--extend-tie-floor "$par_extend_tie_floor"}
  ${par_chain_drop:+-D "$par_chain_drop"}
  ${par_seeded_bases:+-W "$par_seeded_bases"}
  ${par_min_ext_len:+--min-ext-len "$par_min_ext_len"}
  ${par_mate_rescue:+-m "$par_mate_rescue"}
  ${par_rescue_kmer:+--rescue-kmer="$par_rescue_kmer"}
  ${par_rescue_band:+--rescue-band "$par_rescue_band"}
  ${par_rescue_skip:+--rescue-skip}
  ${par_skip_mate_rescue:+-S}
  ${par_skip_pairing:+-P}
  ${par_hic:+--hic}
  ${par_read_type:+-x "$par_read_type"}
  ${par_fast:+--fast}
  ${par_compat:+--compat "$par_compat"}
  ${par_compat_allow_divergent:+--compat-allow-divergent}

  # Scoring options
  ${par_match_score:+-A "$par_match_score"}
  ${par_mismatch_penalty:+-B "$par_mismatch_penalty"}
  ${par_gap_open_penalty:+-O "$par_gap_open_penalty"}
  ${par_gap_extend_penalty:+-E "$par_gap_extend_penalty"}
  ${par_clipping_penalty:+-L "$par_clipping_penalty"}
  ${par_unpaired_penalty:+-U "$par_unpaired_penalty"}

  # Input/Output options
  ${par_smart_pairing:+-p}
  ${par_read_group:+-R "$par_read_group"}
  ${par_header:+-H "$par_header"}
  ${par_ignore_alt:+-j}
  ${par_primary_5prime:+-5}
  ${par_keep_mapq:+-q}
  ${par_batch_size:+-K "$par_batch_size"}
  ${par_chunk_cap:+--chunk-cap "$par_chunk_cap"}
  ${par_cohort_slices:+--cohort-slices "$par_cohort_slices"}
  ${par_cohort_ramp_first:+--cohort-ramp-first "$par_cohort_ramp_first"}
  ${par_cohort_ramp_ratio:+--cohort-ramp-ratio "$par_cohort_ramp_ratio"}
  ${par_verbosity:+-v "$par_verbosity"}
  ${par_min_score:+-T "$par_min_score"}
  ${par_max_hits_xa:+-h "$par_max_hits_xa"}
  ${par_xa_score_fraction:+-z "$par_xa_score_fraction"}
  ${par_output_xb:+-u}
  ${par_output_all:+-a}
  ${par_append_comment:+-C}
  ${par_output_ref_header:+-V}
  ${par_soft_clipping:+-Y}
  ${par_mark_secondary:+-M}
  ${par_insert_size:+-I "$par_insert_size"}

  # Methylation options
  ${par_meth:+--meth="$par_meth"}
  ${par_meth_scoring:+--meth-scoring "$par_meth_scoring"}
  ${par_meth_seed_prune:+--meth-seed-prune="$par_meth_seed_prune"}
  ${par_meth_tags:+--meth-tags="$par_meth_tags"}
  ${par_set_as_failed:+--set-as-failed "$par_set_as_failed"}
  ${par_chimera_qc:+--chimera-qc}

  # Supplementary MAPQ rescoring options
  ${par_supp_rep_hard_cap:+--supp-rep-hard-cap "$par_supp_rep_hard_cap"}
  ${par_proper_pair_from_emitted:+--proper-pair-from-emitted}

  # Seed ordering options
  ${par_seed_order:+--seed-order "$par_seed_order"}

  # Input reader options
  ${par_legacy_reader:+--legacy-reader}

  # Index and input files
  "$par_index"
  "$par_reads1"
  ${par_reads2:+"$par_reads2"}
)

# Run bwa-mem3
bwa-mem3 mem "${cmd_args[@]}"
