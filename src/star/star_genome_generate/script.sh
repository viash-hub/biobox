#!/bin/bash

set -e

## VIASH START
## VIASH END

mkdir -p $par_index

STAR \
    --runMode genomeGenerate \
    --genomeDir $par_index \
    --genomeFastaFiles $par_genome_fasta_files \
    ${meta_cpus:+--runThreadN "${meta_cpus}"} \
    ${par_sjdb_gtf_file:+--sjdbGTFfile "${par_sjdb_gtf_file}"} \
    ${par_sjdbOverhang:+--sjdbOverhang "${par_sjdbOverhang}"} \
    ${par_genome_sa_index_nbases:+--genomeSAindexNbases "${par_genome_sa_index_nbases}"} \
    ${par_sjdb_gtf_chr_prefix:+--sjdbGTFchrPrefix "${par_sjdb_gtf_chr_prefix}"} \
    ${par_sjdb_gtf_feature_exon:+--sjdbGTFfeatureExon "${par_sjdb_gtf_feature_exon}"} \
    ${par_sjdb_gtf_tag_exon_parent_transcript:+--sjdbGTFtag_exon_parent_transcript "${par_sjdb_gtf_tag_exon_parent_transcript}"} \
    ${par_sjdb_gtf_tag_exon_parent_gene:+--sjdbGTFtag_exon_parent_gene "${par_sjdb_gtf_tag_exon_parent_gene}"} \
    ${par_sjdb_gtf_tag_exon_parent_geneName:+--sjdbGTFtag_exon_parent_geneName "${par_sjdb_gtf_tag_exon_parent_geneName}"} \
    ${par_sjdb_gtf_tag_exon_parent_geneType:+--sjdbGTFtag_exon_parent_geneType "${sjdbGTFtag_exon_parent_geneType}"} \
    ${par_limit_genome_generate_ram:+--limitGenomeGenerateRAM "${par_limit_genome_generate_ram}"} \
    ${par_genome_chr_bin_nbits:+--genomeChrBinNbits "${par_genome_chr_bin_nbits}"} \
    ${par_genome_sa_sparse_d:+--genomeSAsparseD "${par_genome_sa_sparse_d}"} \
    ${par_genome_suffix_length_max:+--genomeSuffixLengthMax "${par_genome_suffix_length_max}"} \
    ${par_genome_transform_type:+--genomeTransformType "${par_genome_transform_type}"} \
    ${par_genome_transform_vcf:+--genomeTransformVCF "${par_genome_transform_vCF}"} \
