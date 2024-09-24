#!/bin/bash

set -eo pipefail

unset_if_false=(
    par_sam
    par_error
    par_log2stderr
    par_timeit_header )

for var in "${unset_if_false[@]}"; do
    test_val="${!var}"
    [[ "$test_val" == "false" ]] && unset $var
done

umi_tools prepare-for-rsem \
    ${par_log:+--log "${par_log}"} \
    ${par_tags:+--tags "${par_tags}"} \
    ${par_sam:+--sam} \
    --stdin="${par_input}" \
    ${par_output:+--stdout "${par_output}"} \
    ${par_error:+--error "${par_error}"} \
    ${par_temp_dir:+--temp-dir "${par_temp_dir}"} \
    ${par_log2stderr:+--log2stderr} \
    ${par_verbose:+--verbose "${par_verbose}"} \
    ${par_random_seed:+--random-seed "${par_random_seed}"} \
    ${par_compresslevel:+--compresslevel "${par_compresslevel}"}
    ${par_timeit:+--timeit "${par_timeit}"} \
    ${par_timeit_name:+--timeit-name "${par_timeit_name}"} \
    ${par_timeit_header:+--timeit-header}


