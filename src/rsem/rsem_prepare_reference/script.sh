#!/bin/bash

set -eo pipefail

unset_if_false=( par_gff3_genes_as_transcripts par_polyA par_bowtie par_bowtie2 par_star par_hisat2_hca par_quiet par_prep_pRSEM )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

rsem-prepare-reference \
    ${par_gtf:+--gtf $par_gtf} \
    ${par_gff3:+--gff3 $par_gff3} \
    ${par_gff3_rna_patterns:+--gff3-RNA-patterns $par_gff3_rna_patterns} \
    ${par_gff3_genes_as_transcripts:+--gff3-genes-as-transcripts $par_gff3_genes_as_transcripts} \
    ${par_trusted_sources:+--gtrusted-sourcestf $par_trusted_sources} \
    ${par_transcript_to_gene_map:+--transcript-to-gene-map $par_transcript_to_gene_map} \
    ${par_allele_to_gene_map:+--allele-to-gene-map $par_allele_to_gene_map} \
    ${par_polyA:+--polyA} \
    ${par_polyA_length:+--polyA-length $par_polyA_length} \
    ${par_no_polyA_subset:+--no-polyA-subset $par_no_polyA_subset} \
    ${par_bowtie:+--bowtie} \
    ${par_bowtie_path:+--bowtie-path $par_bowtie_path} \
    ${par_bowtie2:+--bowtie2} \
    ${par_bowtie2_path:+--bowtie2-path $par_bowtie2_path} \
    ${par_star:+--star} \
    ${par_star_path:+--star-path $par_star_path} \
    ${par_star_sjdboverhang:+--star-sjdboverhang $par_star_sjdboverhang} \
    ${par_hisat2_hca:+--hisat2-hca} \
    ${par_hisat2_path:+--hisat2-path $par_hisat2_path} \
    ${par_quiet:+--quiet} \
    ${par_prep_pRSEM:+--prep-pRSEM} \
    ${par_mappability_bigwig_file:+--mappability-bigwig-file $par_mappability_bigwig_file} \
    ${meta_cpus:+--num-threads $meta_cpus} \
    $par_reference_fasta_files \
    $par_reference_name

mkdir -p $par_output
mv $par_reference_name.* $par_output/