#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_inflate" == "false" ]] && unset par_inflate

# Convert a list of file names to multiple -gff arguments
input_files=""
IFS=";" read -ra file_names <<< "$par_gff"
for file in "${file_names[@]}"; do
    input_files+="--gff $file "
done
unset IFS

# take care of --genome (can originally be either a fasta file or an integer)
if [[ -n "$par_genome_size" ]]; then
  genome_arg=$par_genome_size
elif [[ -n "$par_genome_size_fasta" ]]; then
  genome_arg=$par_genome_size_fasta
fi

# run agat_convert_sp_bed2gff.pl
agat_sq_stat_basic.pl \
  $input_files \
  ${genome_arg:+--genome "${genome_arg}"} \
  --output "${par_output}" \
  ${par_inflate:+--inflate} \
  ${par_config:+--config "${par_config}"}
