#!/bin/bash

set -e

## VIASH START
## VIASH END

unset_if_false=(
    par_gencode
    par_features
    par_keep_duplicates
    par_keep_fixed_fasta
    par_sparse
    par_no_clip
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

tmp_dir=$(mktemp -d -p "$meta_temp_dir" "${meta_functionality_name}_XXXXXX")
mkdir -p "$tmp_dir/temp"

if [[ -f "$par_genome" ]] && [[ ! "$par_decoys" ]]; then
    filename="$(basename -- $par_genome)"
    decoys="decoys.txt"
    if [ ${filename##*.} == "gz" ]; then
        grep '^>' <(gunzip -c $par_genome) | cut -d ' ' -f 1 > $decoys
        gentrome="gentrome.fa.gz"
    else
        grep '^>' $par_genome | cut -d ' ' -f 1 > $decoys
        gentrome="gentrome.fa"
    fi
    sed -i.bak -e 's/>//g' $decoys
    cat $par_transcripts $par_genome > $gentrome
else
    gentrome=$par_transcripts
    decoys=$par_decoys
fi

salmon index \
    -t "$gentrome" \
    --tmpdir "$tmp_dir/temp" \
    ${meta_cpus:+--threads "${meta_cpus}"} \
    -i "$par_index" \
    ${par_kmer_len:+-k "${par_kmer_len}"} \
    ${par_gencode:+--gencode} \
    ${par_features:+--features} \
    ${par_keep_duplicates:+--keepDuplicates} \
    ${par_keep_fixed_fasta:+--keepFixedFasta} \
    ${par_filter_size:+-f "${par_filter_size}"} \
    ${par_sparse:+--sparse} \
    ${decoys:+-d "${decoys}"} \
    ${par_no_clip:+--no-clip} \
    ${par_type:+--type "${par_type}"}