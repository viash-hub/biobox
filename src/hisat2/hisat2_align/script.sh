#!/bin/bash

set -euo pipefail

## VIASH START
## VIASH END

# Unset boolean flags that are false
unset_if_false=(
  # Input format
  par_fasta
  par_qseq
  par_phred64
  par_int_quals
  # Paired-end
  par_fr
  par_rf
  par_ff
  par_no_mixed
  par_no_discordant
  par_dovetail
  par_no_contain
  par_no_overlap
  # Alignment
  par_end_to_end
  par_local
  par_very_fast
  par_fast
  par_sensitive
  par_very_sensitive
  par_very_fast_local
  par_fast_local
  par_sensitive_local
  par_very_sensitive_local
  par_ignore_quals
  par_nofw
  par_norc
  par_no_1mm_upfront
  # Spliced alignment
  par_no_spliced_alignment
  par_tmo
  par_dta
  par_dta_cufflinks
  par_avoid_pseudogene
  par_no_templatelen_adjustment
  par_no_temp_splicesite
  par_secondary_alignments
  # Reporting
  par_all
  # SAM
  par_no_unal
  par_no_hd
  par_no_sq
  par_omit_sec_seq
  par_sam_no_qname_trunc
  par_xeq
  par_soft_clipped_unmapped_tlen
  par_add_chrname
  par_remove_chrname
  # Other
  par_new_summary
  par_quiet
  par_reorder
  par_mm
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Build command arguments
cmd_args=(
  ${meta_cpus:+-p "$meta_cpus"}
  -x "$par_index_dir/$par_index_prefix"
  -S "$par_output_sam"
  # Input format
  ${par_fasta:+-f}
  ${par_qseq:+--qseq}
  ${par_phred64:+--phred64}
  ${par_int_quals:+--int-quals}
  ${par_skip:+-s "$par_skip"}
  ${par_qupto:+-u "$par_qupto"}
  ${par_trim5:+--trim5 "$par_trim5"}
  ${par_trim3:+--trim3 "$par_trim3"}
  # Paired-end
  ${par_minins:+-I "$par_minins"}
  ${par_maxins:+-X "$par_maxins"}
  ${par_fr:+--fr}
  ${par_rf:+--rf}
  ${par_ff:+--ff}
  ${par_no_mixed:+--no-mixed}
  ${par_no_discordant:+--no-discordant}
  ${par_dovetail:+--dovetail}
  ${par_no_contain:+--no-contain}
  ${par_no_overlap:+--no-overlap}
  # Alignment
  ${par_end_to_end:+--end-to-end}
  ${par_local:+--local}
  ${par_very_fast:+--very-fast}
  ${par_fast:+--fast}
  ${par_sensitive:+--sensitive}
  ${par_very_sensitive:+--very-sensitive}
  ${par_very_fast_local:+--very-fast-local}
  ${par_fast_local:+--fast-local}
  ${par_sensitive_local:+--sensitive-local}
  ${par_very_sensitive_local:+--very-sensitive-local}
  ${par_n_ceil:+--n-ceil "$par_n_ceil"}
  ${par_dpad:+--dpad "$par_dpad"}
  ${par_gbar:+--gbar "$par_gbar"}
  ${par_ignore_quals:+--ignore-quals}
  ${par_nofw:+--nofw}
  ${par_norc:+--norc}
  ${par_no_1mm_upfront:+--no-1mm-upfront}
  # Spliced alignment
  ${par_no_spliced_alignment:+--no-spliced-alignment}
  ${par_rna_strandness:+--rna-strandness "$par_rna_strandness"}
  ${par_tmo:+--tmo}
  ${par_dta:+--dta}
  ${par_dta_cufflinks:+--dta-cufflinks}
  ${par_avoid_pseudogene:+--avoid-pseudogene}
  ${par_no_templatelen_adjustment:+--no-templatelen-adjustment}
  ${par_known_splicesite_infile:+--known-splicesite-infile "$par_known_splicesite_infile"}
  ${par_novel_splicesite_infile:+--novel-splicesite-infile "$par_novel_splicesite_infile"}
  ${par_novel_splicesite_outfile:+--novel-splicesite-outfile "$par_novel_splicesite_outfile"}
  ${par_no_temp_splicesite:+--no-temp-splicesite}
  ${par_secondary_alignments:+--secondary-alignments}
  ${par_pen_cansplice:+--pen-cansplice "$par_pen_cansplice"}
  ${par_pen_noncansplice:+--pen-noncansplice "$par_pen_noncansplice"}
  ${par_pen_canintronlen:+--pen-canintronlen "$par_pen_canintronlen"}
  ${par_pen_noncanintronlen:+--pen-noncanintronlen "$par_pen_noncanintronlen"}
  ${par_min_intronlen:+--min-intronlen "$par_min_intronlen"}
  ${par_max_intronlen:+--max-intronlen "$par_max_intronlen"}
  # Scoring
  ${par_ma:+--ma "$par_ma"}
  ${par_mp:+--mp "$par_mp"}
  ${par_sp:+--sp "$par_sp"}
  ${par_np:+--np "$par_np"}
  ${par_rdg:+--rdg "$par_rdg"}
  ${par_rfg:+--rfg "$par_rfg"}
  ${par_score_min:+--score-min "$par_score_min"}
  # Reporting
  ${par_k:+-k "$par_k"}
  ${par_all:+-a}
  ${par_max_seeds:+--max-seeds "$par_max_seeds"}
  # SAM
  ${par_no_unal:+--no-unal}
  ${par_no_hd:+--no-hd}
  ${par_no_sq:+--no-sq}
  ${par_rg_id:+--rg-id "$par_rg_id"}
  ${par_omit_sec_seq:+--omit-sec-seq}
  ${par_sam_no_qname_trunc:+--sam-no-qname-trunc}
  ${par_xeq:+--xeq}
  ${par_soft_clipped_unmapped_tlen:+--soft-clipped-unmapped-tlen}
  ${par_add_chrname:+--add-chrname}
  ${par_remove_chrname:+--remove-chrname}
  # Other
  ${par_new_summary:+--new-summary}
  ${par_summary_file:+--summary-file "$par_summary_file"}
  ${par_quiet:+-q}
  ${par_reorder:+--reorder}
  ${par_mm:+--mm}
  ${par_seed:+--seed "$par_seed"}
)

# Handle --rg (multiple values, Viash separator is ;)
if [[ -n "${par_rg:-}" ]]; then
  IFS=';' read -ra _rg_arr <<< "$par_rg"
  for _rg in "${_rg_arr[@]}"; do
    [[ -n "$_rg" ]] && cmd_args+=(--rg "$_rg")
  done
fi

# Single-end vs paired-end mode
if [[ -n "${par_input_r2:-}" ]]; then
  echo "Running in paired-end mode"
  cmd_args+=(-1 "$par_input" -2 "$par_input_r2")
else
  echo "Running in single-end mode"
  cmd_args+=(-U "$par_input")
fi

hisat2 "${cmd_args[@]}"

echo "Alignment complete. Output: $par_output_sam"
