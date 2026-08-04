#!/bin/bash

set -eo pipefail

concat_reads() {
  local out_file="$1"
  shift
  : > "$out_file"
  for f in "$@"; do
    if [ "${f##*.}" == "gz" ]; then
      zcat "$f" >> "$out_file"
    else
      cat "$f" >> "$out_file"
    fi
  done
}

IFS=";" read -ra read_1 <<< "$par_read_1"
IFS=";" read -ra read_2 <<< "$par_read_2"

if [ ${#read_1[@]} -gt 0 ]; then
    concat_reads "$par_fastq_1" "${read_1[@]}"
fi
if [ ${#read_2[@]} -gt 0 ]; then
    concat_reads "$par_fastq_2" "${read_2[@]}"
fi
