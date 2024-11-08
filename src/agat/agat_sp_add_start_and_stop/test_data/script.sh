#!/bin/bash

src=/tmp/agat
dest=src/agat/agat_sp_add_start_and_stop/test_data/

# clone repo
if [ ! -d $src ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT $src
fi

START=300000
END=310000

# subset fasta to 300000 - 310000 range.
seqkit subseq $src/t/scripts_output/in/1.fa -r 300000:310000 > $dest/1.fa


# cp -r /tmp/agat_source/t/scripts_output/in/agat_sp_add_start_and_stop.gff src/agat/agat_sp_add_start_and_stop/test_data/
# cp -r /tmp/agat_source/t/scripts_output/in/1.fa src/agat/agat_sp_add_start_and_stop/test_data/
# cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_add_start_and_stop_1.gff src/agat/agat_sp_add_start_and_stop/test_data/
# cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_add_start_and_stop_2.gff src/agat/agat_sp_add_start_and_stop/test_data/

# subset gff and subtract 300000 from start and end positions
for f in in/agat_sp_add_start_and_stop.gff out/agat_sp_add_start_and_stop_1.gff out/agat_sp_add_start_and_stop_2.gff; do
  src_file=$src/t/scripts_output/$f
  dest_file=$dest/$(basename $f)
  awk -v FS='\t' -v OFS='\t' -v start=$START '
    $1 == "1" && $4 >= start && $5 <= 310000 {
      $4 = $4 - start + 1;
      $5 = $5 - start + 1;
      print;
    }' $src_file > $dest_file
done


# zcat "$ORIG_GTF" | \
#   awk -v FS='\t' -v OFS='\t' "
#     \$1 == \"$CHR\" && \$4 >= $START && \$5 <= $END {
#       \$4 = \$4 - $START + 1;
#       \$5 = \$5 - $START + 1;
#       print;
#     }" > $OUT_DIR/reference_small.gtf