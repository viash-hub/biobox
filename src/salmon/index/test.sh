#!/bin/bash

dir_in="$meta_resources_dir/test_data"

echo ">>> Testing salmon_index"

"$meta_executable" \
  --transcripts "$dir_in/transcriptome.fasta" \
  --index "salmon_index"

echo ">>> Checking whether output exists"
[ ! -d "salmon_index" ] && echo "Salmon index does not exist!" && exit 1
[ -z "$(ls -A 'salmon_index')" ] && echo "Salmon index is empty!" && exit 1

echo "All tests succeeded!"
exit 0
