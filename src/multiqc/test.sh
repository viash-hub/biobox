#!/bin/bash

echo ">>> Testing $meta_functionality_name"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/" \
  --output_report multiqc_report.html \
  --output_data multiqc_data \
  --output_plots multiqc_plots

echo ">>> Check whether output exists"
[ ! -f multiqc_report.html ] && echo "MultiQC report does not exist!" && exit 1
[ ! -d multiqc_data ] && echo "MultiQC data directory does not exist!" && exit 1
[ ! -d multiqc_plots ] && echo "MultiQC plots directory does not exist!" && exit 1

echo "All tests succeeded!"
exit 0