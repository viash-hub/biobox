#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset flags if 'false'
unset_if_false=(
    par_classic
    par_download
    par_file_list
    par_stats
    par_cancer
    par_format_eff
    par_gene_id
    par_hgvs
    par_hgvs_old
    par_hgvs1_letter_aa
    par_hgvs_tr_id
    par_lof
    par_oicr
    par_sequence_ontology
    par_debug
    par_quiet
    par_verbose
    par_canon
    par_interaction
    par_motif
    par_nextprot
    par_only_reg
    par_only_protein
    par_strict
)
for par in ${unset_if_false[@]}; do
    test_val="${!par}" # contains the value of the 'par'
    [[ "$test_val" == "false" ]] && unset $par
done

# Unset flags if 'true'
unset_if_true=(
    par_no_stats
    par_no_downstream
    par_no_intergenic
    par_no_intron
    par_no_upstream
    par_no_utr
    par_no_hgvs
    par_no_lof
    par_no_shift_hgvs
    par_no_download
    par_no_log
    par_no_tag
    par_no_genome
    par_no_expand_iub
    par_no_interaction
    par_no_motif
    par_no_nextprot
)
for par in ${unset_if_true[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "true" ]] && unset $par
done

# Run SnpEff
snpEff \
    ${par_chr:+-chr "$par_chr"} \
    ${par_classic:+-classic} \
    ${par_csv_stats:+-csvStats "$par_csv_stats"} \
    ${par_download:+-download} \
    ${par_input_format:+-i "$par_input_format"} \
    ${par_file_list:+-fileList} \
    ${par_output_format:+-o "$par_output_format"} \
    ${par_stats:+-stats} \
    ${par_no_stats:+-noStats} \
    ${par_fi:+-fi "$par_fi"} \
    ${par_no_downstream:+-no-downstream} \
    ${par_no_intergenic:+-no-intergenic} \
    ${par_no_intron:+-no-intron} \
    ${par_no_upstream:+-no-upstream} \
    ${par_no_utr:+-no-utr} \
    ${par_no:+-no "$par_no"} \
    ${par_cancer:+-cancer} \
    ${par_cancer_samples:+-cancerSamples "$par_cancer_samples]"} \
    ${par_fastaprot:+-fastaProt "$par_fastaprot]"} \
    ${par_format_eff:+-formatEff} \
    ${par_gene_id:+-geneId} \
    ${par_hgvs:+-hgvs} \
    ${par_hgvs_old:+-hgvsOld} \
    ${par_hgvs1_letter_aa:+-hgvs1LetterAa} \
    ${par_hgvs_tr_id:+-hgvsTrId} \
    ${par_lof:+-lof} \
    ${par_no_hgvs:+-noHgvs} \
    ${par_no_lof:+-noLof} \
    ${par_no_shift_hgvs:+-noShiftHgvs} \
    ${par_oicr:+-oicr} \
    ${par_sequence_ontology:+-sequenceOntology} \
    ${par_config:+-config "$par_config"} \
    ${par_config_option:+-configOption "$par_config_option"} \
    ${par_debug:+-debug} \
    ${par_data_dir:+-dataDir "$par_data_dir"} \
    ${par_no_download:+-nodownload} \
    ${par_no_log:+-noLog} \
    ${par_quiet:+-quiet} \
    ${par_verbose:+-verbose} \
    ${par_canon:+-canon} \
    ${par_canon_list:+-canonList "$par_canon_list"} \
    ${par_tag:+-tag "$par_tag"} \
    ${par_no_tag:+-notag} \
    ${par_interaction:+-interaction} \
    ${par_interval:+-interval "$par_interval"} \
    ${par_max_tsl:+-maxTSL "$par_max_tsl"} \
    ${par_motif:+-motif} \
    ${par_nextprot:+-nextProt} \
    ${par_no_genome:+-noGenome} \
    ${par_no_expand_iub:+-noExpandIUB} \
    ${par_no_interaction:+-noInteraction} \
    ${par_no_motif:+-noMotif} \
    ${par_no_nextprot:+-noNextProt} \
    ${par_only_reg:+-onlyReg} \
    ${par_only_protein:+-onlyProtein} \
    ${par_only_tr:+-onlyTr "$par_onlyTr"} \
    ${par_reg:+-reg "$par_reg"} \
    ${par_ss:+-ss "$par_ss"} \
    ${par_splice_region_exon_size:+-spliceRegionExonSize "$par_splice_region_exon_size"} \
    ${par_splice_region_intron_min:+-spliceRegionIntronMin "$par_splice_region_intron_min"} \
    ${par_splice_region_intron_max:+-spliceRegionIntronMax "$par_splice_region_intron_max"} \
    ${par_strict:+-strict} \
    ${par_ud:+-ud "$par_ud"} \
    "$par_genome_version" \
    "$par_input" \
    > "$par_output"

# Path of the output file (par_output)
absolute_path=$(realpath "$par_output")
directory_path=$(dirname "$absolute_path")

# Move the automatically generated outputs to the pwd
if [ -z "$par_no_stats" ]; then
    mv -f snpEff_genes.txt snpEff_summary.html "$directory_path"
fi