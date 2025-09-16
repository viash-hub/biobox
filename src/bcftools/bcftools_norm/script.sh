#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
unset_if_false=(
  par_atomize
  par_remove_duplicates_flag
  par_force
  par_no_version
  par_do_not_normalize
  par_strict_filter
)

for par in ${unset_if_false[@]}; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Build command array
cmd_args=(
  bcftools norm
  ${par_atomize:+--atomize}
  ${par_atom_overlaps:+--atom-overlaps "$par_atom_overlaps"}
  ${par_check_ref:+--check-ref "$par_check_ref"}
  ${par_remove_duplicates_flag:+--remove-duplicates}
  ${par_rm_dup:+--rm-dup "$par_rm_dup"}
  ${par_exclude:+--exclude "$par_exclude"}
  ${par_fasta_ref:+--fasta-ref "$par_fasta_ref"}
  ${par_force:+--force}
  ${par_gff_annot:+--gff-annot "$par_gff_annot"}
  ${par_include:+--include "$par_include"}
  ${par_keep_sum:+--keep-sum "$par_keep_sum"}
  ${par_multiallelics:+--multiallelics "$par_multiallelics"}
  ${par_multi_overlaps:+--multi-overlaps "$par_multi_overlaps"}
  ${par_no_version:+--no-version}
  ${par_do_not_normalize:+--do-not-normalize}
  ${par_old_rec_tag:+--old-rec-tag "$par_old_rec_tag"}
  ${par_output_type:+--output-type "$par_output_type"}
  ${par_regions:+--regions "$par_regions"}
  ${par_regions_file:+--regions-file "$par_regions_file"}
  ${par_regions_overlap:+--regions-overlap "$par_regions_overlap"}
  ${par_strict_filter:+--strict-filter}
  ${par_sort:+--sort "$par_sort"}
  ${par_targets:+--targets "$par_targets"}
  ${par_targets_file:+--targets-file "$par_targets_file"}
  ${par_targets_overlap:+--targets-overlap "$par_targets_overlap"}
  ${meta_cpus:+--threads "$meta_cpus"}
  ${par_verbosity:+--verbosity "$par_verbosity"}
  ${par_site_win:+--site-win "$par_site_win"}
  ${par_write_index:+--write-index="$par_write_index"}
  ${par_output:+--output "$par_output"}
  "$par_input"
)

# Execute command
"${cmd_args[@]}"
