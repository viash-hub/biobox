#!/bin/bash

set -eo pipefail

filename="$(basename -- "$par_input")"

if [ ${filename##*.} == "gz" ]; then
    gunzip -c "$par_input" > "$par_output"
else
    cat "$par_input" > "$par_output"
fi
