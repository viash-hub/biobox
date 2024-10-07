#!/bin/bash

echo ">>> Testing $meta_executable"

test_dir="${meta_resources_dir}/test_data"

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/rsem.tar.gz
# gunzip -k rsem.tar.gz
# tar -xf rsem.tar
# mv $test_dir/rsem $meta_resources_dir

echo "> Prepare test data"

cat > reads_R1.fastq <<'EOF'
@SEQ_ID1
ACGCTGCCTCATAAGCCTCACACAT
+
IIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
ACCCGCAAGATTAGGCTCCGTACAC
+
!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > reads_R2.fastq <<'EOF'
@SEQ_ID1
ATGTGTGAGGCTTATGAGGCAGCGT
+
IIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
GTGTACGGAGCCTAATCTTGCAGGG
+
!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > genome.fasta <<'EOF'
>chr1
TGGCATGAGCCAACGAACGCTGCCTCATAAGCCTCACACATCCGCGCCTATGTTGTGACTCTCTGTGAGCGTTCGTGGG
GCTCGTCACCACTATGGTTGGCCGGTTAGTAGTGTGACTCCTGGTTTTCTGGAGCTTCTTTAAACCGTAGTCCAGTCAA
TGCGAATGGCACTTCACGACGGACTGTCCTTAGGTGTGAGGCTTATGAGGCACTCAGGGGA
EOF

cat > genes.gtf <<'EOF'
chr1	example_source	gene	0	50	.	+	.	gene_id "gene1"; transcript_id "transcript1";
chr1	example_source	exon	20	40	.	+	.	gene_id "gene1"; transcript_id "transcript1";
chr1	example_source	gene	100	219	.	+	.	gene_id "gene2"; transcript_id "transcript2";
chr1	example_source	exon	191	210	.	+	.	gene_id "gene2"; transcript_id "transcript2";
EOF

cat > ref.cnt <<'EOF'
1 0 0 1
0 0 0
0 3
0	1
Inf	0
EOF

cat > ref.genes.results <<'EOF'
gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM
gene1	transcript1	21.00	21.00	0.00	0.00	0.00
gene2	transcript2	20.00	20.00	0.00	0.00	0.00
EOF

cat > ref.isoforms.results <<'EOF'
transcript_id	gene_id	length	effective_length	expected_count	TPM	FPKM	IsoPct
transcript1	gene1	21	21.00	0.00	0.00	0.00	0.00
transcript2	gene2	20	20.00	0.00	0.00	0.00	0.00
EOF


echo "> Generate index"

rsem-prepare-reference \
  --gtf "genes.gtf" \
  "genome.fasta" \
  "index"

mkdir index
mv index.* index/
  
STAR \
  ${meta_cpus:+--runThreadN $meta_cpus} \
  --runMode genomeGenerate \
  --genomeDir "index/" \
  --genomeFastaFiles "genome.fasta" \
  --sjdbGTFfile "genes.gtf" \
  --genomeSAindexNbases 2
  
#########################################################################################

echo ">>> Test 1: Paired-end reads using STAR to align reads"
"$meta_executable" \
	--star \
	--paired \
	--input "reads_R1.fastq;reads_R2.fastq" \
	--index index \
	--id test \
	--seed 1 \
	--quiet

echo ">>> Checking whether output exists"
[ ! -f "test.genes.results" ] && echo "Gene level expression counts file does not exist!" && exit 1
[ ! -s "test.genes.results" ] && echo "Gene level expression counts file is empty!" && exit 1
[ ! -f "test.isoforms.results" ] && echo "Transcript level expression counts file does not exist!" && exit 1
[ ! -s "test.isoforms.results" ] && echo "Transcript level expression counts file is empty!" && exit 1
[ ! -d "test.stat" ] && echo "Stats file does not exist!" && exit 1

echo ">>> Check wheter output is correct"
diff ref.genes.results test.genes.results || { echo "Gene level expression counts file is incorrect!"; exit 1; }
diff ref.isoforms.results test.isoforms.results || { echo "Transcript level expression counts file is incorrect!"; exit 1; }
diff ref.cnt test.stat/test.cnt || { echo "Stats file is incorrect!"; exit 1; }

#####################################################################################################

echo "All tests succeeded!"
exit 0
