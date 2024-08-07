#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/scripts_output/in/1.gff src/agat/agat_sp_compare_two_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/in/agat_sp_compare_two_annotations/file1.gff src/agat/agat_sp_compare_two_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/in/agat_sp_compare_two_annotations/file2.gff src/agat/agat_sp_compare_two_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_compare_two_annotations_1.txt src/agat/agat_sp_compare_two_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_compare_two_annotations_2.txt src/agat/agat_sp_compare_two_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_compare_two_annotations_3.txt src/agat/agat_sp_compare_two_annotations/test_data