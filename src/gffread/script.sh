#!/bin/bash

## VIASH START
## VIASH END


gffread \
    --input "$par_input" \
    ${par_chr_mapping:+-m "$par_chr_mapping"} \
    ${par_seq_info:+-s "$par_seq_info"} \
    ${par_outfile:+-o "$par_outfile"} \
    ${par_force_exons:+--force_exons} \
    ${par_gene2exon:+--gene2exon} \
    ${par_t_adopt:+--t-adopt} \
    ${par_decode:+-D} \
    ${par_merge_exons:+-Z} \
    ${par_genome:+-g} \
    ${par_junctions:+-j} \
    ${par_spliced_exons:+-w} \
    ${par_w_add:+--w_add "$par_w_add"} \
    ${par_w_nocds:+--w_nocds} \
    ${par_spliced_cds:+-x} \
    ${par_tr_cds:+-y} \
    ${par_w_coords:+-W} \
    ${par_stop_dot:+-S} \
    ${par_id_version:+-L} \
    ${par_trackname:+-t "$par_trackname"} \
    ${par_gtf_output:+-T} \
    ${par_bed:+--bed} \
    ${par_tlf:+--tlf} \
    ${par_table:+--table} \
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
    ${par_ignore_locus:+--ignore_locus} \
    ${par_description:+-A} \
    ${par_sort_alpha:+--sort_alpha} \
    ${par_sort_by:+--sort_by "$par_sort_by"} \
    ${par_keep_attrs:+-F} \
    ${par_keep_exon_attrs:+--keep_exon_attrs} \
    ${par_no_exon_attrs:+-G} \
    ${par_attrs:+--attrs "$par_attrs"} \
    ${par_keep_genes:+--keep_genes} \
    ${par_keep_comments:+--keep_comments} \
    ${par_process_other:+-O} \
    ${par_rm_stop_codons:+-V} \
    ${par_adj_cds_start:+-H} \
    ${par_opposite_strand:+-B} \
    ${par_coding_status:+-P} \
    ${par_add_hasCDS:+--add_hasCDS} \
    ${par_adj_stop:+--adj_stop} \
    ${par_rm_noncanon:+-N} \
    ${par_complete_cds:+-J} \
    ${par_no_pseudo:+--no_pseudo} \
    ${par_in_bed:+--in_bed} \
    ${par_in_tlf:+--in_tlf} \
    ${par_stream:+--stream} \
    ${par_merge:+-M} \
    ${par_dupinfo:+-d "$par_dupinfo"} \
    ${par_cluster_only:+--cluster_only} \
    ${par_rm_redundant:+-K} \
    ${par_no_boundary:+-Q} \
    ${par_no_overlap:+-Y} \

