#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_coding" == "false" ]] && unset par_coding
[[ "$par_strict_range" == "false" ]] && unset par_strict_range
[[ "$par_no_single_exon" == "false" ]] && unset par_no_single_exon
[[ "$par_nc" == "false" ]] && unset par_nc
[[ "$par_ignore_locus" == "false" ]] && unset par_ignore_locus
[[ "$par_description" == "false" ]] && unset par_description
[[ "$par_sort_alpha" == "false" ]] && unset par_sort_alpha
[[ "$par_keep_attrs" == "false" ]] && unset par_keep_attrs
[[ "$par_keep_exon_attrs" == "false" ]] && unset par_keep_exon_attrs
[[ "$par_keep_comments" == "false" ]] && unset par_keep_comments
[[ "$par_process_other" == "false" ]] && unset par_process_other
[[ "$par_rm_stop_codons" == "false" ]] && unset par_rm_stop_codons
[[ "$par_adj_cds_start" == "false" ]] && unset par_adj_cds_start
[[ "$par_opposite_strand" == "false" ]] && unset par_opposite_strand
[[ "$par_coding_status" == "false" ]] && unset par_coding_status
[[ "$par_add_hasCDS" == "false" ]] && unset par_add_hasCDS
[[ "$par_adj_stop" == "false" ]] && unset par_adj_stop
[[ "$par_rm_noncanon" == "false" ]] && unset par_rm_noncanon
[[ "$par_complete_cds" == "false" ]] && unset par_complete_cds
[[ "$par_no_pseudo" == "false" ]] && unset par_no_pseudo
[[ "$par_in_bed" == "false" ]] && unset par_in_bed
[[ "$par_in_tlf" == "false" ]] && unset par_in_tlf
[[ "$par_stream" == "false" ]] && unset par_stream
[[ "$par_merge" == "false" ]] && unset par_merge
[[ "$par_rm_redundant" == "false" ]] && unset par_rm_redundant
[[ "$par_no_boundary" == "false" ]] && unset par_no_boundary
[[ "$par_no_overlap" == "false" ]] && unset par_no_overlap
[[ "$par_force_exons" == "false" ]] && unset par_force_exons
[[ "$par_gene2exon" == "false" ]] && unset par_gene2exon
[[ "$par_t_adopt" == "false" ]] && unset par_t_adopt
[[ "$par_decode" == "false" ]] && unset par_decode
[[ "$par_merge_exons" == "false" ]] && unset par_merge_exons
[[ "$par_junctions" == "false" ]] && unset par_junctions
[[ "$par_w_nocds" == "false" ]] && unset par_w_nocds
[[ "$par_tr_cds" == "false" ]] && unset par_tr_cds
[[ "$par_w_coords" == "false" ]] && unset par_w_coords
[[ "$par_stop_dot" == "false" ]] && unset par_stop_dot
[[ "$par_id_version" == "false" ]] && unset par_id_version
[[ "$par_gtf_output" == "false" ]] && unset par_gtf_output
[[ "$par_bed" == "false" ]] && unset par_bed
[[ "$par_tlf" == "false" ]] && unset par_tlf
[[ "$par_expose_dups" == "false" ]] && unset par_expose_dups
[[ "$par_cluster_only" == "false" ]] && unset par_cluster_only


gffread \
    "$par_input" \
    ${par_chr_mapping:+-m "$par_chr_mapping"} \
    ${par_seq_info:+-s "$par_seq_info"} \
    ${par_outfile:+-o "$par_outfile"} \
    ${par_force_exons:+--force-exons} \
    ${par_gene2exon:+--gene2exon} \
    ${par_t_adopt:+--t-adopt} \
    ${par_decode:+-D} \
    ${par_merge_exons:+-Z} \
    ${par_genome:+-g "$par_genome"} \
    ${par_junctions:+-j} \
    ${par_spliced_exons:+-w "$par_spliced_exon"} \
    ${par_w_add:+--w-add "$par_w_add"} \
    ${par_w_nocds:+--w-nocds} \
    ${par_spliced_cds:+-x "$par_spliced_cds"} \
    ${par_tr_cds:+-y "$par_tr_cds"} \
    ${par_w_coords:+-W} \
    ${par_stop_dot:+-S} \
    ${par_id_version:+-L} \
    ${par_trackname:+-t "$par_trackname"} \
    ${par_gtf_output:+-T} \
    ${par_bed:+--bed} \
    ${par_tlf:+--tlf} \
    ${par_table:+--table "$par_table"} \
    ${par_expose_dups:+-E} \
    ${par_ids:+--ids "$par_ids"} \
    ${par_nids:+--nids "$par_nids"} \
    ${par_maxintron:+-i "$par_maxintron"} \
    ${par_minlen:+-l "$par_minlen"} \
    ${par_range:+-r "$par_range"} \
    ${par_strict_range:+-R} \
    ${par_jmatch:+--jmatch "$par_jmatch"} \
    ${par_no_single_exon:+-U} \
    ${par_coding:+-C} \
    ${par_nc:+--nc} \
    ${par_ignore_locus:+--ignore-locus} \
    ${par_description:+-A} \
    ${par_sort_alpha:+--sort-alpha} \
    ${par_sort_by:+--sort-by "$par_sort_by"} \
    ${par_keep_attrs:+-F} \
    ${par_keep_exon_attrs:+--keep-exon-attrs} \
    ${par_no_exon_attrs:+-G} \
    ${par_attrs:+--attrs "$par_attrs"} \
    ${par_keep_genes:+--keep-genes} \
    ${par_keep_comments:+--keep-comments} \
    ${par_process_other:+-O} \
    ${par_rm_stop_codons:+-V} \
    ${par_adj_cds_start:+-H} \
    ${par_opposite_strand:+-B} \
    ${par_coding_status:+-P} \
    ${par_add_hasCDS:+--add-hasCDS} \
    ${par_adj_stop:+--adj-stop} \
    ${par_rm_noncanon:+-N} \
    ${par_complete_cds:+-J} \
    ${par_no_pseudo:+--no-pseudo} \
    ${par_in_bed:+--in-bed} \
    ${par_in_tlf:+--in-tlf} \
    ${par_stream:+--stream} \
    ${par_merge:+-M} \
    ${par_dupinfo:+-d "$par_dupinfo"} \
    ${par_cluster_only:+--cluster-only} \
    ${par_rm_redundant:+-K} \
    ${par_no_boundary:+-Q} \
    ${par_no_overlap:+-Y} 

