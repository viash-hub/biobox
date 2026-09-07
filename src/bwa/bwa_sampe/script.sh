#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_preload_index" == "false" ]] && unset par_preload_index
[[ "$par_disable_smith_waterman" == "false" ]] && unset par_disable_smith_waterman
[[ "$par_disable_insert_size_estimate" == "false" ]] && unset par_disable_insert_size_estimate

# Resolve the index base name from the index directory. BWA takes a prefix that
# its five index files share, which is not itself a path that can be mounted or
# staged, so the directory holding them is passed in instead.
index_dir="${par_index%/}"

if [[ ! -d "$index_dir" ]]; then
    echo "Error: --index must be a directory containing the BWA index files," >&2
    echo "but '$index_dir' is not a directory." >&2
    exit 1
fi

mapfile -d '' -t bwt_files < <(find -L "$index_dir" -maxdepth 1 -name '*.bwt' -type f -print0)

if [[ "${#bwt_files[@]}" -eq 0 ]]; then
    echo "Error: No BWA index found in '$index_dir': expected a '.bwt' file." >&2
    exit 1
elif [[ "${#bwt_files[@]}" -gt 1 ]]; then
    echo "Error: Multiple BWA indices found in '$index_dir': ${bwt_files[*]}." >&2
    echo "The index directory must contain exactly one index." >&2
    exit 1
fi

index_prefix="${bwt_files[0]%.bwt}"

# Build the command
cmd_args=(
    # Pairing options
    ${par_max_insert_size:+-a "$par_max_insert_size"}
    ${par_max_occ_one_end:+-o "$par_max_occ_one_end"}
    ${par_max_hits_paired:+-n "$par_max_hits_paired"}
    ${par_max_hits_discordant:+-N "$par_max_hits_discordant"}
    ${par_chimeric_rate:+-c "$par_chimeric_rate"}
    
    # Output options
    ${par_output:+-f "$par_output"}
    ${par_read_group:+-r "$par_read_group"}
    
    # Algorithm options
    ${par_preload_index:+-P}
    ${par_disable_smith_waterman:+-s}
    ${par_disable_insert_size_estimate:+-A}
    
    # Required arguments: index, SAI files, FASTQ files
    "$index_prefix"
    "$par_sai1"
    "$par_sai2"
    "$par_reads1"
    "$par_reads2"
)

# Run bwa sampe
bwa sampe "${cmd_args[@]}"
