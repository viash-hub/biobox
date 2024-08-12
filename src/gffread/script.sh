#!/bin/bash

## VIASH START
## VIASH END

# unset flags
unset_if_false=(
    par_coding
    par_strict_range
    par_no_single_exon
    par_no_exon_attrs
    par_nc
    par_ignore_locus
    par_description
    par_sort_alpha
    par_keep_genes
    par_keep_attrs
    par_keep_exon_attrs
    par_keep_comments
    par_process_other
    par_rm_stop_codons
    par_adj_cds_start
    par_opposite_strand
    par_coding_status
    par_add_hasCDS
    par_adj_stop
    par_rm_noncanon
    par_complete_cds
    par_no_pseudo
    par_in_bed
    par_in_tlf
    par_stream
    par_merge
    par_rm_redundant
    par_no_boundary
    par_no_overlap
    par_force_exons
    par_gene2exon
    par_t_adopt
    par_decode
    par_merge_exons
    par_junctions
    par_w_nocds
    par_tr_cds
    par_w_coords
    par_stop_dot
    par_id_version
    par_gtf_output
    par_bed
    par_tlf
    par_expose_dups
    par_cluster_only
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# if par_table is not empty, replace ";" with ","
par_table=$(echo "$par_table" | tr ';' ',')

$(which gffread) \
    "$par_input" \
    ${par_chr_mapping:+-m "$par_chr_mapping"} \
    ${par_seq_info:+-s "$par_seq_info"} \
    -o "$par_outfile" \
    ${par_force_exons:+--force-exons} \
    ${par_gene2exon:+--gene2exon} \
    ${par_t_adopt:+--t-adopt} \
    ${par_decode:+-D} \
    ${par_merge_exons:+-Z} \
    ${par_genome:+-g "$par_genome"} \
    ${par_junctions:+-j} \
    ${par_spliced_exons:+-w "$par_spliced_exons"} \
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

