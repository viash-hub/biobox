#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_skip_duplicate_marking" == "false" ]] && unset par_skip_duplicate_marking
[[ "$par_extra_information" == "false" ]] && unset par_extra_information
[[ "$par_fill_gaps" == "false" ]] && unset par_fill_gaps

# replace ';' with ','
par_interesting_contigs=$(echo $par_interesting_contigs | tr ';' ',')
par_viral_contigs=$(echo $par_viral_contigs | tr ';' ',')
par_disable_filters=$(echo $par_disable_filters | tr ';' ',')

# run arriba
arriba \
  -x "$par_bam" \
  -a "$par_genome" \
  -g "$par_gene_annotation" \
  -o "$par_fusions" \
  ${par_known_fusions:+-k "${par_known_fusions}"} \
  ${par_blacklist:+-b "${par_blacklist}"} \
  ${par_structural_variants:+-d "${par_structural_variants}"} \
  ${par_tags:+-t "${par_tags}"} \
  ${par_protein_domains:+-p "${par_protein_domains}"} \
  ${par_fusions_discarded:+-O "${par_fusions_discarded}"} \
  ${par_max_genomic_breakpoint_distance:+-D "${par_max_genomic_breakpoint_distance}"} \
  ${par_strandedness:+-s "${par_strandedness}"} \
  ${par_interesting_contigs:+-i "${par_interesting_contigs}"} \
  ${par_viral_contigs:+-v "${par_viral_contigs}"} \
  ${par_disable_filters:+-f "${par_disable_filters}"} \
  ${par_max_e_value:+-E "${par_max_e_value}"} \
  ${par_min_supporting_reads:+-S "${par_min_supporting_reads}"} \
  ${par_max_mismappers:+-m "${par_max_mismappers}"} \
  ${par_max_homolog_identity:+-L "${par_max_homolog_identity}"} \
  ${par_homopolymer_length:+-H "${par_homopolymer_length}"} \
  ${par_read_through_distance:+-R "${par_read_through_distance}"} \
  ${par_min_anchor_length:+-A "${par_min_anchor_length}"} \
  ${par_many_spliced_events:+-M "${par_many_spliced_events}"} \
  ${par_max_kmer_content:+-K "${par_max_kmer_content}"} \
  ${par_max_mismatch_pvalue:+-V "${par_max_mismatch_pvalue}"} \
  ${par_fragment_length:+-F "${par_fragment_length}"} \
  ${par_max_reads:+-U "${par_max_reads}"} \
  ${par_quantile:+-Q "${par_quantile}"} \
  ${par_exonic_fraction:+-e "${par_exonic_fraction}"} \
  ${par_top_n:+-T "${par_top_n}"} \
  ${par_covered_fraction:+-C "${par_covered_fraction}"} \
  ${par_max_itd_length:+-l "${par_max_itd_length}"} \
  ${par_min_itd_allele_fraction:+-z "${par_min_itd_allele_fraction}"} \
  ${par_min_itd_supporting_reads:+-Z "${par_min_itd_supporting_reads}"} \
  ${par_skip_duplicate_marking:+-u} \
  ${par_extra_information:+-X} \
  ${par_fill_gaps:+-I}
