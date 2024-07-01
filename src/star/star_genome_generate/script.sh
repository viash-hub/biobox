#!/bin/bash

set -e

## VIASH START
## VIASH END

mkdir -p $par_index

STAR \
    --runMode genomeGenerate \
    --genomeDir $par_index \
    --genomeFastaFiles $par_genomeFastaFiles \
    ${meta_cpus:+--runThreadN "${meta_cpus}"} \
    ${par_sjdbGTFfile:+--sjdbGTFfile "${par_sjdbGTFfile}"} \
    ${par_sjdbOverhang:+--sjdbOverhang "${par_sjdbOverhang}"} \
    ${par_genomeSAindexNbases:+--genomeSAindexNbases "${par_genomeSAindexNbases}"} \
    ${par_sjdbGTFchrPrefix:+--sjdbGTFchrPrefix "${par_sjdbGTFchrPrefix}"} \
    ${par_sjdbGTFfeatureExon:+--sjdbGTFfeatureExon "${par_sjdbGTFfeatureExon}"} \
    ${par_sjdbGTFtagExonParentTranscript:+--sjdbGTFtagExonParentTranscript "${par_sjdbGTFtagExonParentTranscript}"} \
    ${par_sjdbGTFtagExonParentGene:+--sjdbGTFtagExonParentGene "${par_sjdbGTFtagExonParentGene}"} \
    ${par_sjdbGTFtagExonParentGeneName:+--sjdbGTFtagExonParentGeneName "${par_sjdbGTFtagExonParentGeneName}"} \
    ${par_sjdbGTFtagExonParentGeneType:+--sjdbGTFtagExonParentGeneType "${sjdbGTFtagExonParentGeneType}"} \
    ${par_limitGenomeGenerateRAM:+--limitGenomeGenerateRAM "${par_limitGenomeGenerateRAM}"} \
    ${par_genomeChrBinNbits:+--genomeChrBinNbits "${par_genomeChrBinNbits}"} \
    ${par_genomeSAsparseD:+--genomeSAsparseD "${par_genomeSAsparseD}"} \
    ${par_genomeSuffixLengthMax:+--genomeSuffixLengthMax "${par_genomeSuffixLengthMax}"} \
    ${par_genomeTransformType:+--genomeTransformType "${par_genomeTransformType}"} \
    ${par_genomeTransformVCF:+--genomeTransformVCF "${par_genomeTransformVCF}"} \
