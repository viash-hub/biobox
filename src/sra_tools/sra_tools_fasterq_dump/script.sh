#!/usr/bin/env bash

unset_if_false=( 
    par_details
    par_progress
    par_split_spot
    par_split_files
    par_split_3
    par_concatenate_reads
    skip_technical
    include_technical
    exclude_technical
    par_fasta
    par_fasta_unsorted
    par_fasta_reference_table
    par_fasta_concat_all
    par_internal_ref
    par_external_ref
    par_only_unaligned
    par_only_aligned
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

if [ -z "$par_accession" ] && [ -z "$par_prefetch_directory" ]; then
    echo "Either 'accesssion' or 'prefetch_directory' must be specified."
    exit 1
fi

if [ ! -z "$par_accession" ] && [ ! -z "$par_prefetch_directory" ]; then
    echo "'accesssion' or 'prefetch_directory' are mutually exclusive arguments."
    exit 1
fi
input=${par_accession-$par_prefetch_directory}



fasterq-dump \
    ${par_details:+--details} \
    ${par_progress:+--progress} \
    ${par_split_spot:+--split-spot} \
    ${par_split_files:+--split-files} \
    ${par_concatenate_reads:+--concatenate-reads} \
    ${par_split_3:+--split-3} \
    ${par_skip_technical:+--skip-technical} \
    ${par_include_technical:+--include-technical} \
    ${par_exclude_technical:+--exclude-technical} \
    ${par_minimal_read_length:+--min-read-len $par_minimal_read_length} \
    ${par_table:+--table $par_table} \
    ${par_bases:+--bases $par_bases} \
    ${par_fasta+--fasta} \
    ${par_fasta_unsorted+--fasta-unsorted} \
    ${par_fasta_reference_table+--fasta-ref-tbl} \
    ${par_fasta_concat_all+--fasta-concat-all} \
    ${par_internal_ref+--internal-ref} \
    ${par_external_ref+--external-ref} \
    ${par_seq_defline:+--seq-defline ${par_seq_defline@Q}} \
    ${par_qual_defline:+--qual-defline ${par_seq_defline@Q}} \
    ${par_only_unaligned+--only-unaligned} \
    ${par_only_aligned+--only-aligned} \
    ${par_log_level+--log-level $par_log_level} \
    ${meta_memory_mb:+--memory "${meta_memory_mb}M"} \
    --force \
    "$input"
    