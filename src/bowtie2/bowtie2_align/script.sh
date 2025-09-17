#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_fastq" == "false" ]] && unset par_fastq
[[ "$par_tab5" == "false" ]] && unset par_tab5
[[ "$par_tab6" == "false" ]] && unset par_tab6
[[ "$par_qseq" == "false" ]] && unset par_qseq
[[ "$par_fasta" == "false" ]] && unset par_fasta
[[ "$par_raw" == "false" ]] && unset par_raw
[[ "$par_cmdline" == "false" ]] && unset par_cmdline
[[ "$par_phred33" == "false" ]] && unset par_phred33
[[ "$par_phred64" == "false" ]] && unset par_phred64
[[ "$par_int_quals" == "false" ]] && unset par_int_quals
[[ "$par_very_fast" == "false" ]] && unset par_very_fast
[[ "$par_fast" == "false" ]] && unset par_fast
[[ "$par_sensitive" == "false" ]] && unset par_sensitive
[[ "$par_very_sensitive" == "false" ]] && unset par_very_sensitive
[[ "$par_very_fast_local" == "false" ]] && unset par_very_fast_local
[[ "$par_fast_local" == "false" ]] && unset par_fast_local
[[ "$par_sensitive_local" == "false" ]] && unset par_sensitive_local
[[ "$par_very_sensitive_local" == "false" ]] && unset par_very_sensitive_local
[[ "$par_ignore_quals" == "false" ]] && unset par_ignore_quals
[[ "$par_nofw" == "false" ]] && unset par_nofw
[[ "$par_norc" == "false" ]] && unset par_norc
[[ "$par_no_1mm_upfront" == "false" ]] && unset par_no_1mm_upfront
[[ "$par_end_to_end" == "false" ]] && unset par_end_to_end
[[ "$par_local" == "false" ]] && unset par_local
[[ "$par_all" == "false" ]] && unset par_all
[[ "$par_fr" == "false" ]] && unset par_fr
[[ "$par_rf" == "false" ]] && unset par_rf
[[ "$par_ff" == "false" ]] && unset par_ff
[[ "$par_no_mixed" == "false" ]] && unset par_no_mixed
[[ "$par_no_discordant" == "false" ]] && unset par_no_discordant
[[ "$par_dovetail" == "false" ]] && unset par_dovetail
[[ "$par_no_contain" == "false" ]] && unset par_no_contain
[[ "$par_no_overlap" == "false" ]] && unset par_no_overlap
[[ "$par_time" == "false" ]] && unset par_time
[[ "$par_quiet" == "false" ]] && unset par_quiet
[[ "$par_met_stderr" == "false" ]] && unset par_met_stderr
[[ "$par_no_unal" == "false" ]] && unset par_no_unal
[[ "$par_no_head" == "false" ]] && unset par_no_head
[[ "$par_no_sq" == "false" ]] && unset par_no_sq
[[ "$par_omit_sec_seq" == "false" ]] && unset par_omit_sec_seq
[[ "$par_sam_no_qname_trunc" == "false" ]] && unset par_sam_no_qname_trunc
[[ "$par_xeq" == "false" ]] && unset par_xeq
[[ "$par_soft_clipped_unmapped_tlen" == "false" ]] && unset par_soft_clipped_unmapped_tlen
[[ "$par_sam_append_comment" == "false" ]] && unset par_sam_append_comment
[[ "$par_align_paired_reads" == "false" ]] && unset par_align_paired_reads
[[ "$par_preserve_tags" == "false" ]] && unset par_preserve_tags
[[ "$par_reorder" == "false" ]] && unset par_reorder
[[ "$par_mm" == "false" ]] && unset par_mm
[[ "$par_qc_filter" == "false" ]] && unset par_qc_filter
[[ "$par_non_deterministic" == "false" ]] && unset par_non_deterministic

# Validate input arguments
if [[ -z "$par_index" ]]; then
  echo "Error: --index is required" >&2
  exit 1
fi

# Validate that at least one input type is specified
if [[ -z "$par_mate1" && -z "$par_mate2" && -z "$par_unpaired" && -z "$par_interleaved" && -z "$par_bam_input" ]]; then
  echo "Error: At least one input type must be specified (--mate1/--mate2, --unpaired, --interleaved, or --bam_input)" >&2
  exit 1
fi

# Validate paired-end input
if [[ -n "$par_mate1" && -z "$par_mate2" ]] || [[ -z "$par_mate1" && -n "$par_mate2" ]]; then
  echo "Error: Both --mate1 and --mate2 must be specified for paired-end reads" >&2
  exit 1
fi

# Build the command arguments
cmd_args=(
    -x "$par_index"
    ${par_mate1:+-1 "$(IFS=','; echo "${par_mate1[*]}")"}
    ${par_mate2:+-2 "$(IFS=','; echo "${par_mate2[*]}")"}
    ${par_unpaired:+-U "$(IFS=','; echo "${par_unpaired[*]}")"}
    ${par_interleaved:+--interleaved "$(IFS=','; echo "${par_interleaved[*]}")"}
    ${par_bam_input:+-b "$(IFS=','; echo "${par_bam_input[*]}")"}
    -S "$par_output"
    ${par_fastq:+-q}
    ${par_tab5:+--tab5}
    ${par_tab6:+--tab6}
    ${par_qseq:+--qseq}
    ${par_fasta:+-f}
    ${par_raw:+-r}
    ${par_cmdline:+-c}
    ${par_skip:+-s "$par_skip"}
    ${par_upto:+-u "$par_upto"}
    ${par_trim5:+-5 "$par_trim5"}
    ${par_trim3:+-3 "$par_trim3"}
    ${par_trim_to:+--trim-to "$par_trim_to"}
    ${par_continuous_fasta:+-F "$par_continuous_fasta"}
    ${par_phred33:+--phred33}
    ${par_phred64:+--phred64}
    ${par_int_quals:+--int-quals}
    ${par_very_fast:+--very-fast}
    ${par_fast:+--fast}
    ${par_sensitive:+--sensitive}
    ${par_very_sensitive:+--very-sensitive}
    ${par_very_fast_local:+--very-fast-local}
    ${par_fast_local:+--fast-local}
    ${par_sensitive_local:+--sensitive-local}
    ${par_very_sensitive_local:+--very-sensitive-local}
    ${par_N:+-N "$par_N"}
    ${par_L:+-L "$par_L"}
    ${par_i:+-i "$par_i"}
    ${par_n_ceil:+--n-ceil "$par_n_ceil"}
    ${par_dpad:+--dpad "$par_dpad"}
    ${par_gbar:+--gbar "$par_gbar"}
    ${par_ignore_quals:+--ignore-quals}
    ${par_nofw:+--nofw}
    ${par_norc:+--norc}
    ${par_no_1mm_upfront:+--no-1mm-upfront}
    ${par_end_to_end:+--end-to-end}
    ${par_local:+--local}
    ${par_ma:+--ma "$par_ma"}
    ${par_mp:+--mp "$par_mp"}
    ${par_np:+--np "$par_np"}
    ${par_rdg:+--rdg "$par_rdg"}
    ${par_rfg:+--rfg "$par_rfg"}
    ${par_score_min:+--score-min "$par_score_min"}
    ${par_k:+-k "$par_k"}
    ${par_all:+-a}
    ${par_D:+-D "$par_D"}
    ${par_R:+-R "$par_R"}
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
    ${par_time:+-t}
    ${par_un:+--un "$par_un"}
    ${par_al:+--al "$par_al"}
    ${par_un_conc:+--un-conc "$par_un_conc"}
    ${par_al_conc:+--al-conc "$par_al_conc"}
    ${par_quiet:+--quiet}
    ${par_met_file:+--met-file "$par_met_file"}
    ${par_met_stderr:+--met-stderr}
    ${par_met:+--met "$par_met"}
    ${par_no_unal:+--no-unal}
    ${par_no_head:+--no-head}
    ${par_no_sq:+--no-sq}
    ${par_rg_id:+--rg-id "$par_rg_id"}
    ${par_rg:+--rg "$par_rg"}
    ${par_omit_sec_seq:+--omit-sec-seq}
    ${par_sam_no_qname_trunc:+--sam-no-qname-trunc}
    ${par_xeq:+--xeq}
    ${par_soft_clipped_unmapped_tlen:+--soft-clipped-unmapped-tlen}
    ${par_sam_append_comment:+--sam-append-comment}
    ${par_sam_opt_config:+--sam-opt-config "$par_sam_opt_config"}
    ${par_align_paired_reads:+--align-paired-reads}
    ${par_preserve_tags:+--preserve-tags}
    ${meta_cpus:+-p "$meta_cpus"}
    ${par_reorder:+--reorder}
    ${par_mm:+--mm}
    ${par_qc_filter:+--qc-filter}
    ${par_seed:+--seed "$par_seed"}
    ${par_non_deterministic:+--non-deterministic}
)

# Run bowtie2
bowtie2 "${cmd_args[@]}"
