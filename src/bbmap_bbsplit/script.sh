#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT 

if [ ! -d "$par_built_bbsplit_index" ]; then
    other_refs=()
    while IFS="," read -r name path 
    do
        other_refs+=("ref_$name=$path")
    done < "$par_bbsplit_fasta_list"
fi

if $par_only_build_index; then
    if [ -f "$par_primary_ref" ] && [ ${#other_refs[@]} -gt 0 ]; then
        bbsplit.sh \
            ref_primary="$par_primary_ref" "${other_refs[@]}" \
            path=$par_bbsplit_index \
            threads=${meta_cpus:-1}
    else
        echo "ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files."
    fi
else
    IFS="," read -ra input <<< "$par_input"
    tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
    index_files=''
    if [ -d "$par_built_bbsplit_index" ]; then
        index_files="path=$par_built_bbsplit_index"
    elif [ -f "$par_primary_ref" ] && [ ${#other_refs[@]} -gt 0 ]; then
        index_files="ref_primary=$par_primary_ref ${other_refs[@]}"
    else
        echo "ERROR: Please either specify a BBSplit index as input or a primary fasta file along with names and paths to non-primary fasta files."
    fi
    if $par_paired; then
        bbsplit.sh \
            $index_files \
            threads=${meta_cpus:-1} \
            in=${input[0]} \
            in2=${input[1]} \
            basename=${tmpdir}/%_#.fastq \
            refstats=bbsplit_stats.txt
        read1=$(find $tmpdir/ -iname primary_1*)
        read2=$(find $tmpdir/ -iname primary_2*)
        cp $read1 $par_fastq_1
        cp $read2 $par_fastq_2
    else
        bbsplit.sh \
            $index_files \
            threads=${meta_cpus:-1} \
            in=${input[0]} \
            basename=${tmpdir}/%.fastq \
            refstats=bbsplit_stats.txt
        read1=$(find $tmpdir/ -iname primary*)
        cp $read1 $par_fastq_1
    fi
fi

exit 0