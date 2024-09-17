#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset flags if 'false'
unset_if_false=(
    par_classic
    par_download
    par_fileList
    par_stats
    par_cancer
    par_formatEff
    par_geneId
    par_hgvs
    par_hgvsOld
    par_hgvs1LetterAa
    par_hgvsTrId
    par_lof
    par_oicr
    par_sequenceOntology
    par_debug
    par_quiet
    par_verbose
    par_canon
    par_interaction
    par_motif
    par_nextProt
    par_onlyReg
    par_onlyProtein
    par_strict
)
for par in ${unset_if_false[@]}; do
    test_val="${!par}" # contains the value of the 'par'
    [[ "$test_val" == "false" ]] && unset $par
done

# Unset flags if 'true'
unset_if_true=(
    par_noStats
    par_no_downstream
    par_no_intergenic
    par_no_intron
    par_no_upstream
    par_no_utr
    par_noHgvs
    par_noLof
    par_noShiftHgvs
    par_nodownload
    par_noLog
    par_notag
    par_noGenome
    par_noExpandIUB
    par_noInteraction
    par_noMotif
    par_noNextProt
)
for par in ${unset_if_true[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "true" ]] && unset $par
done

# Run SnpEff
snpEff \
    ${par_chr:+-chr "$par_chr"} \
    ${par_classic:+-classic} \
    ${par_csvStats:+-csvStats "$par_csvStats"} \
    ${par_download:+-download} \
    ${par_input_format:+-i "$par_input_format"} \
    ${par_fileList:+-fileList} \
    ${par_output_format:+-o "$par_output_format"} \
    ${par_stats:+-stats} \
    ${par_noStats:+-noStats} \
    ${par_fi:+-fi "$par_fi"} \
    ${par_no_downstream:+-no-downstream} \
    ${par_no_intergenic:+-no-intergenic} \
    ${par_no_intron:+-no-intron} \
    ${par_no_upstream:+-no-upstream} \
    ${par_no_utr:+-no-utr} \
    ${par_no:+-no "$par_no"} \
    ${par_cancer:+-cancer} \
    ${par_cancerSamples:+-cancerSamples "$par_cancerSamples]"} \
    ${par_fastaProt:+-fastaProt "$par_fastaProt]"} \
    ${par_formatEff:+-formatEff} \
    ${par_geneId:+-geneId} \
    ${par_hgvs:+-hgvs} \
    ${par_hgvsOld:+-hgvsOld} \
    ${par_hgvs1LetterAa:+-hgvs1LetterAa} \
    ${par_hgvsTrId:+-hgvsTrId} \
    ${par_lof:+-lof} \
    ${par_noHgvs:+-noHgvs} \
    ${par_noLof:+-noLof} \
    ${par_noShiftHgvs:+-noShiftHgvs} \
    ${par_oicr:+-oicr} \
    ${par_sequenceOntology:+-sequenceOntology} \
    ${par_config:+-config "$par_config"} \
    ${par_configOption:+-configOption "$par_configOption"} \
    ${par_debug:+-debug} \
    ${par_dataDir:+-dataDir "$par_dataDir"} \
    ${par_nodownload:+-nodownload} \
    ${par_noLog:+-noLog} \
    ${par_quiet:+-quiet} \
    ${par_verbose:+-verbose} \
    ${par_canon:+-canon} \
    ${par_canonList:+-canonList "$par_canonList"} \
    ${par_tag:+-tag "$par_tag"} \
    ${par_notag:+-notag} \
    ${par_interaction:+-interaction} \
    ${par_interval:+-interval "$par_interval"} \
    ${par_maxTSLl:+-maxTSL "$par_maxTSL"} \
    ${par_motif:+-motif} \
    ${par_nextProt:+-nextProt} \
    ${par_noGenome:+-noGenome} \
    ${par_noExpandIUB:+-noExpandIUB} \
    ${par_noInteraction:+-noInteraction} \
    ${par_noMotif:+-noMotif} \
    ${par_noNextProt:+-noNextProt} \
    ${par_onlyReg:+-onlyReg} \
    ${par_onlyProtein:+-onlyProtein} \
    ${par_onlyTr:+-onlyTr "$par_onlyTr"} \
    ${par_reg:+-reg "$par_reg"} \
    ${par_ss:+-ss "$par_ss"} \
    ${par_spliceRegionExonSize:+-spliceRegionExonSize "$par_spliceRegionExonSize"} \
    ${par_spliceRegionIntronMin:+-spliceRegionIntronMin "$par_spliceRegionIntronMin"} \
    ${par_spliceRegionIntronMax:+-spliceRegionIntronMax "$par_spliceRegionIntronMax"} \
    ${par_strict:+-strict} \
    ${par_ud:+-ud "$par_ud"} \
    "$par_genome_version" \
    "$par_input" \
    > "$par_output"

# Path of the output file (par_output)
absolute_path=$(realpath "$par_output")
directory_path=$(dirname "$absolute_path")

# Move the automatically generated outputs to the pwd
if [ -z "$par_noStats" ]; then
    mv -f snpEff_genes.txt snpEff_summary.html "$directory_path"
fi
# mv -f snpEff_genes.txt "$directory_path"