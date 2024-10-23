#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

unset_if_false=( par_fastx par_sq par_fastx par_print_all_reads par_paired_in par_paired_out
                 par_F par_R par_verbose par_de_novo par_otu_map par_full_search par_out2
                 par_sout par_sam par_paired )


for var in "${unset_if_false[@]}"; do
    if [ "${!var}" == "false" ]; then
        unset $var
    fi
done

reads=()
IFS=";" read -ra input <<< "$par_input"
if [ "${#input[@]}" -eq 2 ]; then
    reads="--reads ${input[0]} --reads ${input[1]}"
    # set paired to true in case it's not
    par_paired=true
else
    reads="--reads ${input[0]}"
    par_paired=false
fi

refs=()

# check if references are input normally or through a manifest file
if [[ ! -z "$par_ribo_database_manifest" ]]; then
    while IFS= read -r path || [[ -n $path ]]; do
        refs=$refs" --ref $path"
    done < $par_ribo_database_manifest

elif [[ ! -z "$par_ref" ]]; then
    IFS=";" read -ra ref <<< "$par_ref"
    # check if length is 2 and par_paired is set to true
    if [[ "${#ref[@]}" -eq 2 && "$par_paired" == "true" ]]; then
        refs="--ref ${ref[0]} --ref ${ref[1]}"
    # check if length is 1 and par_paired is set to false
    elif [[ "${#ref[@]}" -eq 1 && "$par_paired" == "false" ]]; then
            refs="--ref $par_ref"      
    else # if one reference provided but paired is set to true:
        echo "Two reference fasta files are required for paired-end reads"
            exit 1
    fi
else 
    echo "No reference fasta file(s) provided"
    exit 1
fi


sortmerna \
    $refs \
    $reads \
    ${par_output:+--aligned "${par_output}"} \
    ${par_fastx:+--fastx} \
    ${par_other:+--other "${par_other}"} \
    ${par_kvdb:+--kvdb "${par_kvdb}"} \
    ${par_idx_dir:+--idx-dir "${par_idx_dir}"} \
    ${par_readb:+--readb "${par_readb}"} \
    ${par_sam:+--sam} \
    ${par_sq:+--sq} \
    ${par_blast:+--blast "${par_blast}"} \
    ${par_num_alignments:+--num_alignments "${par_num_alignments}"} \
    ${par_min_lis:+--min_lis "${par_min_lis}"} \
    ${par_print_all_reads:+--print_all_reads} \
    ${par_paired_in:+--paired_in} \
    ${par_paired_out:+--paired_out} \
    ${par_out2:+--out2} \
    ${par_sout:+--sout} \
    ${par_zip_out:+--zip-out "${par_zip_out}"} \
    ${par_match:+--match "${par_match}"} \
    ${par_mismatch:+--mismatch "${par_mismatch}"} \
    ${par_gap_open:+--gap_open "${par_gap_open}"} \
    ${par_gap_ext:+--gap_ext "${par_gap_ext}"} \
    ${par_N:+-N "${par_N}"} \
    ${par_a:+-a "${par_a}"} \
    ${par_e:+-e "${par_e}"} \
    ${par_F:+-F} \
    ${par_R:+-R} \
    ${par_num_alignment:+--num_alignment "${par_num_alignment}"} \
    ${par_best:+--best "${par_best}"} \
    ${par_verbose:+--verbose} \
    ${par_id:+--id "${par_id}"} \
    ${par_coverage:+--coverage "${par_coverage}"} \
    ${par_de_novo:+--de_novo} \
    ${par_otu_map:+--otu_map} \
    ${par_num_seed:+--num_seed "${par_num_seed}"} \
    ${par_passes:+--passes "${par_passes}"} \
    ${par_edge:+--edge "${par_edge}"} \
    ${par_full_search:+--full_search} \
    ${par_index:+--index "${par_index}"} \
    ${par_L:+-L $par_L} \
    ${par_interval:+--interval "${par_interval}"} \
    ${par_max_pos:+--max_pos "${par_max_pos}"}


if [ ! -z $par_log ]; then
    mv "${par_output}.log" $par_log
fi

exit 0

