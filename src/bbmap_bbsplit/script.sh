#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT 

unset_if_false=( par_paired par_only_build_index par_interleaved par_untrim par_nzo)

for var in "${unset_if_false[@]}"; do
    if [ -z "${!var}" ]; then
        unset $var
    fi
done

if [ ! -d "$par_build" ]; then
    IFS=";" read -ra ref_files <<< "$par_ref"
    primary_ref="${ref_files[0]}"
    refs=()
    for file in "${ref_files[@]:1}"
    do
        name=$(basename "$file" | sed 's/\.[^.]*$//')
        refs+=("ref_$name=$file")
    done
fi

if $par_only_build_index; then
    if [ ${#refs[@]} -gt 1 ]; then
        bbsplit.sh \
            --ref_primary="$primary_ref" \
            "${refs[@]}" \
            path=$par_build
    else
        echo "ERROR: Please specify at least two reference fasta files."
    fi
else
    IFS=";" read -ra input <<< "$par_input"
    tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
    index_files=''
    if [ -d "$par_build" ]; then
        index_files="path=$par_build"
    elif [ ${#refs[@]} -gt 0 ]; then
        index_files="--ref_primary=$primary_ref ${refs[*]}"
    else
        echo "ERROR: Please either specify a BBSplit index as input or at least two reference fasta files."
    fi

    extra_args=""
    if [ -n "$par_refstats" ]; then extra_args+=" --refstats $par_refstats"; fi
    if [ -n "$par_ambiguous" ]; then extra_args+=" --ambiguous $par_ambiguous"; fi
    if [ -n "$par_ambiguous2" ]; then extra_args+=" --ambiguous2 $par_ambiguous2"; fi
    if [ -n "$par_minratio" ]; then extra_args+=" --minratio $par_minratio"; fi
    if [ -n "$par_minhits" ]; then extra_args+=" --minhits $par_minhits"; fi
    if [ -n "$par_maxindel" ]; then extra_args+=" --maxindel $par_maxindel"; fi
    if [ -n "$par_qin" ]; then extra_args+=" --qin $par_qin"; fi
    if [ -n "$par_qtrim" ]; then extra_args+=" --qtrim $par_qtrim"; fi
    if [ "$par_interleaved" = true ]; then extra_args+=" --interleaved"; fi
    if [ "$par_untrim" = true ]; then extra_args+=" --untrim"; fi
    if [ "$par_nzo" = true ]; then extra_args+=" --nzo"; fi

    if [ -n "$par_bbmap_args" ]; then extra_args+=" $par_bbmap_args"; fi

    
    if $par_paired; then
        bbsplit.sh \
            $index_files \
            in=${input[0]} \
            in2=${input[1]} \
            basename=${tmpdir}/%_#.fastq \
            $extra_args
        read1=$(find $tmpdir/ -iname primary_1*)
        read2=$(find $tmpdir/ -iname primary_2*)
        cp $read1 $par_fastq_1
        cp $read2 $par_fastq_2
    else
        bbsplit.sh \
            $index_files \
            in=${input[0]} \
            basename=${tmpdir}/%.fastq \
            $extra_args
        read1=$(find $tmpdir/ -iname primary*)
        cp $read1 $par_fastq_1
    fi
fi

exit 0
