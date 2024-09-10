#!/bin/bash

set -eo pipefail

python3 "$meta_resources_dir/prepare-for-rsem.py" \
    ${par_log:+--log "${par_log}"} \
    ${par_tags:+--tags "${par_tags}"} \
    ${par_sam:+--sam} \
    --stdin="${par_input}" \
    ${par_output:+--stdout "${par_output}"}
