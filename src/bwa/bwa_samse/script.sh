#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

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
    # Options
    ${par_max_occ:+-n "$par_max_occ"}
    ${par_output:+-f "$par_output"}
    ${par_read_group:+-r "$par_read_group"}
    
    # Required arguments: index, SAI file, FASTQ file
    "$index_prefix"
    "$par_sai"
    "$par_reads"
)

# Run bwa samse
bwa samse "${cmd_args[@]}"
