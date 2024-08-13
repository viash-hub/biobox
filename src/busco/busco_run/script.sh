#!/bin/bash

## VIASH START
## VIASH END

unset_if_false=(
    par_tar
    par_force
    par_quiet
    par_restart
    par_auto_lineage
    par_auto_lineage_euk
    par_auto_lineage_prok
    par_augustus
    par_long
    par_scaffold_composition
    par_miniprot
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

tmp_dir=$(mktemp -d -p "$meta_temp_dir" busco_XXXXXXXXX)
prefix=$(openssl rand -hex 8)

busco \
    --in "$par_input" \
    --mode "$par_mode" \
    --out "$prefix" \
    --out_path "$tmp_dir" \
    --opt-out-run-stats \
    ${meta_cpus:+--cpu "${meta_cpus}"} \
    ${par_lineage_dataset:+--lineage_dataset "$par_lineage_dataset"} \
    ${par_augustus:+--augustus} \
    ${par_augustus_parameters:+--augustus_parameters "$par_augustus_parameters"} \
    ${par_augustus_species:+--augustus_species "$par_augustus_species"} \
    ${par_auto_lineage:+--auto-lineage} \
    ${par_auto_lineage_euk:+--auto-lineage-euk} \
    ${par_auto_lineage_prok:+--auto-lineage-prok} \
    ${par_contig_break:+--contig_break $par_contig_break} \
    ${par_datasets_version:+--datasets_version "$par_datasets_version"} \
    ${par_e_value:+--evalue "$par_e_value"} \
    ${par_force:+--force} \
    ${par_limit:+--limit "$par_limit"} \
    ${par_long:+--long} \
    ${par_metaeuk:+--metaeuk} \
    ${par_metaeuk_parameters:+--metaeuk_parameters "$par_metaeuk_parameters"} \
    ${par_metaeuk_rerun_parameters:+--metaeuk_rerun_parameters "$par_metaeuk_rerun_parameters"} \
    ${par_miniprot:+--miniprot} \
    ${par_quiet:+--quiet} \
    ${par_restart:+--restart} \
    ${par_scaffold_composition:+--scaffold_composition} \
    ${par_tar:+--tar} \


out_dir=$(find "$tmp_dir/$prefix" -maxdepth 1 -name 'run_*')

if [[ -n "$par_short_summary_json" ]]; then
    cp "$out_dir/short_summary.json" "$par_short_summary_json"
fi
if [[ -n "$par_short_summary_txt" ]]; then
    cp "$out_dir/short_summary.txt" "$par_short_summary_txt"
fi
if [[ -n "$par_full_table" ]]; then
    cp "$out_dir/full_table.tsv" "$par_full_table"
fi
if [[ -n "$par_missing_busco_list" ]]; then
    cp "$out_dir/missing_busco_list.tsv" "$par_missing_busco_list"
fi
if [[ -n "$par_output_dir" ]]; then
    if [[ -d "$par_output_dir" ]]; then
        rm -r "$par_output_dir"
    fi
    cp -r -L "$out_dir" "$par_output_dir"
fi

