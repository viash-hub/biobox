#!/bin/bash

echo ">>> Test $meta_functionality_name"

echo "> Prepare test data"

cat > reads_R1.fastq <<'EOF'
@SEQ_ID1
ACAGGGTTTCACCATGTTGGCCAGG
+
IIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
TCCCAGGTAACAAACCAACCAACTT
+
!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > reads_R2.fastq <<'EOF'
@SEQ_ID1
TACCATTACCCTACCATCCACCATG
+
IIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
CACTCGGCTGCATGCTTAGTGCACT
+
!!!!!!!!!!!!!!!!!!!!!!!!!
EOF

cat > genome.fasta <<'EOF'
>I
AGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTGGTCTTGATCTCCTGACCTCAGGTGATCCATCCGCCT
TGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCACCTGGCCTGGTTTCGAACTCTTGACCTCAGGTGGTCTG
CCCATCTTGACCTTCCAAAGTGCTGGAGCTACAGGCATGAGCCACTGCACCTGGTGCTTTTGGTAAAAGCAACCTGGAAT
CAAATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTT
TAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGAC
EOF

cat > human.fa <<'EOF'
>human
AGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTGGTCTTGATCTCCTGACCTCAGGTGATCCATCCGCCT
TGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCACCTGGCCTGGTTTCGAACTCTTGACCTCAGGTGGTCTG
CCCATCTTGACCTTCCAAAGTGCTGGAGCTACAGGCATGAGCCACTGCACCTGGTGCTTTTGGTAAAAGCAACCTGGAAT
EOF

cat > sarscov2.fa <<'EOF'
>sarscov2
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAA
AATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGG
ACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTT
EOF

####################################################################################################

echo ">>> Building BBSplit index"
"${meta_executable}" \
  --ref "genome.fasta;human.fa;sarscov2.fa" \
  --only_build_index \
  --build "BBSplit_index" 

echo ">>> Check whether output exists"
[ ! -d "BBSplit_index" ] && echo "BBSplit index does not exist!" && exit 1
[ -z "$(ls -A 'BBSplit_index')" ] && echo "BBSplit index is empty!" && exit 1

####################################################################################################


echo ">>> Testing with single-end reads and primary/non-primary FASTA files"
"${meta_executable}" \
  --input "reads_R1.fastq" \
  --ref "genome.fasta;human.fa;sarscov2.fa" \
  --fastq_1 "filtered_reads_R1.fastq"

echo ">>> Check whether output exists"
[ ! -f "filtered_reads_R1.fastq" ] && echo "Filtered reads file does not exist!" && exit 1
[ ! -s "filtered_reads_R1.fastq" ] && echo "Filtered reads file is empty!" && exit 1

echo ">>> Check whether output is correct"
grep -q "ACAGGGTTTCACCATGTTGGCCAGG" filtered_reads_R1.fastq || { echo "Filtered reads file does not contain expected sequence!"; exit 1; }

rm filtered_reads_R1.fastq

####################################################################################################

echo ">>> Testing with paired-end reads and primary/non-primary FASTA files"
"${meta_executable}" \
  --paired \
  --input "reads_R1.fastq;reads_R2.fastq" \
  --ref "genome.fasta;human.fa;sarscov2.fa" \
  --fastq_1 "filtered_reads_R1.fastq" \
  --fastq_2 "filtered_reads_R2.fastq"

echo ">>> Check whether output exists"
[ ! -f "filtered_reads_R1.fastq" ] && echo "Filtered read 1 file does not exist!" && exit 1
[ ! -s "filtered_reads_R1.fastq" ] && echo "Filtered read 1 file is empty!" && exit 1
[ ! -f "filtered_reads_R2.fastq" ] && echo "Filtered read 2 file does not exist!" && exit 1
[ ! -s "filtered_reads_R2.fastq" ] && echo "Filtered read 2 file is empty!" && exit 1

echo ">>> Check whether output is correct"
grep -q "ACAGGGTTTCACCATGTTGGCCAGG" filtered_reads_R1.fastq || { echo "Filtered read 1 file does not contain expected sequence!"; exit 1; }
grep -q "TACCATTACCCTACCATCCACCATG" filtered_reads_R2.fastq || { echo "Filtered read 2 file does not contain expected sequence!"; exit 1; }

rm filtered_reads_R1.fastq filtered_reads_R2.fastq

####################################################################################################

echo ">>> Testing with single-end reads and BBSplit index"
"${meta_executable}" \
  --input "reads_R1.fastq" \
  --build "BBSplit_index" \
  --fastq_1 "filtered_reads_R1.fastq"

echo ">>> Check whether output exists"
[ ! -f "filtered_reads_R1.fastq" ] && echo "Filtered reads file does not exist!" && exit 1
[ ! -s "filtered_reads_R1.fastq" ] && echo "Filtered reads file is empty!" && exit 1

echo ">>> Check whether output is correct"
grep -q "ACAGGGTTTCACCATGTTGGCCAGG" filtered_reads_R1.fastq || { echo "Filtered reads file does not contain expected sequence!"; exit 1; }

rm filtered_reads_R1.fastq

####################################################################################################

echo ">>> Testing with paired-end reads and BBSplit index"
"${meta_executable}" \
  --paired \
  --input "reads_R1.fastq;reads_R2.fastq" \
  --build "BBSplit_index" \
  --fastq_1 "filtered_reads_R1.fastq" \
  --fastq_2 "filtered_reads_R2.fastq"

echo ">>> Check whether output exists"
[ ! -f "filtered_reads_R1.fastq" ] && echo "Filtered read 1 file does not exist!" && exit 1
[ ! -s "filtered_reads_R1.fastq" ] && echo "Filtered read 1 file is empty!" && exit 1
[ ! -f "filtered_reads_R2.fastq" ] && echo "Filtered read 2 file does not exist!" && exit 1
[ ! -s "filtered_reads_R2.fastq" ] && echo "Filtered read 2 file is empty!" && exit 1


echo ">>> Check whether output is correct"
grep -q "ACAGGGTTTCACCATGTTGGCCAGG" filtered_reads_R1.fastq || { echo "Filtered read 1 file does not contain expected sequence!"; exit 1; }
grep -q "TACCATTACCCTACCATCCACCATG" filtered_reads_R2.fastq || { echo "Filtered read 2 file does not contain expected sequence!"; exit 1; }

rm filtered_reads_R1.fastq filtered_reads_R2.fastq

echo "All tests succeeded!"
exit 0