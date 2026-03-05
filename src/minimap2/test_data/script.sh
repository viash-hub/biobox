# test_data_creation

## minimap2_index

### --input "$meta_resources_dir/test_data/test.fasta"
cat <<EOF > "$meta_resources_dir/test_data/test.fasta"
>chr1
ACGTACGTACGT
EOF

### check.mmi TODO, no cat creation/read possible!
#cat <<EOF > "$meta_resources_dir/test_data/check.mmi"
#MMI
#chr1
#    222
#EOF


## minimap_align

### --reference "$TEST_DATA/ref.fasta" \
cat <<EOF > "$TEST_DATA/ref.fasta"
>seq1
ACTGATCGATCGATCGATCGATCGATCGATCGATCGATCGACTATCGATCGATCGATCGA
EOF

### --query "$TEST_DATA/reads.fastq" \
cat <<EOF > "$TEST_DATA/reads.fastq"
@read1
ACTGATCGATCGATCGATCGATCAAAAGATCGATCGATCGACTATCGATCGATCGATCGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

### check-paf
cat <<EOF > "$TEST_DATA/check.paf"
read1	60	5	60	+	seq1	60	9	60	45	55	6	tp:A:P	cm:i:6	s1:i:44	s2:i:0	dv:f:0	rl:i:0
EOF

### check-bam
cat <<EOF > "$TEST_DATA/check.bam"
read1	60	5	60	+	seq1	60	9	60	45	55	6	tp:A:P	cm:i:6	s1:i:44	s2:i:0	dv:f:0	rl:i:0
EOF

