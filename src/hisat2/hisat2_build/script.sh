#!/bin/bash

set -euo pipefail

## VIASH START
## VIASH END

# Unset boolean flags that are false
unset_if_false=(
  par_large_index
  par_noauto
  par_nodc
  par_noref
  par_justref
  par_quiet
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

mkdir -p "$par_index_dir"

cmd_args=(
  ${meta_cpus:+-p "$meta_cpus"}
  ${par_large_index:+--large-index}
  ${par_noauto:+--noauto}
  ${par_nodc:+--nodc}
  ${par_noref:+--noref}
  ${par_justref:+--justref}
  ${par_bmax:+--bmax "$par_bmax"}
  ${par_bmaxdivn:+--bmaxdivn "$par_bmaxdivn"}
  ${par_dcv:+--dcv "$par_dcv"}
  ${par_offrate:+--offrate "$par_offrate"}
  ${par_ftabchars:+--ftabchars "$par_ftabchars"}
  ${par_localoffrate:+--localoffrate "$par_localoffrate"}
  ${par_localftabchars:+--localftabchars "$par_localftabchars"}
  ${par_ss:+--ss "$par_ss"}
  ${par_exon:+--exon "$par_exon"}
  ${par_snp:+--snp "$par_snp"}
  ${par_haplotype:+--haplotype "$par_haplotype"}
  ${par_seed:+--seed "$par_seed"}
  ${par_quiet:+-q}
  "$par_reference"
  "$par_index_dir/$par_index_prefix"
)

hisat2-build "${cmd_args[@]}"

echo "hisat2-build complete. Index written to: $par_index_dir/$par_index_prefix"
