#!/bin/bash

set -eo pipefail

if $par_star; then 
    STAR --runMode genomeGenerate \
        --genomeDir $par_rsem \
        --genomeFastaFiles $par_fasta \
        --sjdbGTFfile $par_gtf \
        ${meta_cpus:+--runThreadN $meta_cpus} \
        ${meta_memory:+--limitGenomeGenerateRAM $meta_memory}
    rsem-prepare-reference \
        --gtf $par_gtf \
        --star \
        ${meta_cpus:+--num-threads $meta_cpus} \
        $par_fasta \
        $par_rsem/genome
    cp "$par_rsem/genome.transcripts.fa" $par_transcript_fasta
else 
    rsem-prepare-reference \
    --gtf $par_gtf \
    ${meta_cpus:+--num-threads $meta_cpus} \
    $par_fasta \
    "$par_rsem/genome"
    cp "$par_rsem/genome.transcripts.fa" $par_transcript_fasta
fi

# Version
text="${meta_functionality_name}:
    rsem: $(rsem-calculate-expression --version | sed -e 's/Current version: RSEM v//g')
    star: $(STAR --version | sed -e 's/STAR_//g')"
if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi     