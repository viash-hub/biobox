#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo ">>> Run salmon_index"
"$meta_executable" \
  --transcripts $dir_in/transcriptome.fasta \
  --index index

echo ">>> Checking whether output exists"
[ ! -d "index" ] && echo "'index' does not exist!" && exit 1
[ -z "$(ls -A 'index')" ] && echo "'index' is empty!" && exit 1

echo "All tests succeeded!"
exit 0
