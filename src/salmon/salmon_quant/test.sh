#!/bin/bash

set -e

echo "==============================================================================="
echo "> Prepare test data"

dir_in="test_data"
mkdir -p "$dir_in"

cat > "$dir_in/transcriptome.fasta" <<'EOF'
>contig1
AGCTCCAGATTCGCTCAGGCCCTTGATCATCAGTCGTCGTCGTCTTCGATTTGCCAGAGG
AGTTTAGATGAAGAATGTCAAGGATGTTCCTCCCTGCCCTCCCATCTAGCCAAGAACATT
TCCAAGAAGATAAAACTGTCACTGAGACAGGTCTGGATGCGCCCTAGGGGCAAATAGAGA
>contig2
AGGCCTTTACCACATTGCTGCTGGCTATAGGAAGTCCCAGGTACTAGCCTGAAACAGCTG
ATATTTGGGGCTGTCACAGACAATATGGCCACCCCTTGGTCTTTATGCATGAAGATTATG
TAAAGGTTTTTATTAAAAAATATATATATATATATAAATGATCTAGATTATTTTCCTCTT
TCTGAAGTACTTTCTTAAAAAAATAAAATTAAATGTTTATAGTATTCCCGGT
EOF

cat > "$dir_in/a_1.fq" <<'EOF'
@SEQ_ID1
AGAATGTCAAGGATGTTCCTCC
+
IIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
ACCCGCAAGATTAGGCTCCGTA
+
!!!!!!!!!!!!!!!!!!!!!!
@SEQ_ID3
CTCAGGCCCTTGATCATCAGTC
+
IIIIIIIIIIIIIIIIIIIIII
EOF

cat > "$dir_in/a_2.fq" <<'EOF'
@SEQ_ID1
GGAGGAACATCCTTGACATTCT
+
IIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
GTGTACGGAGCCTAATCTTGCA
+
!!!!!!!!!!!!!!!!!!!!!!
@SEQ_ID3
GACTGATGATCAAGGGCCTGAG
+
IIIIIIIIIIIIIIIIIIIIII
EOF

cat > "$dir_in/b_1.fq" <<'EOF'
@SEQ_ID1
CTTTACCACATTGCTGCTGGCT
+
IIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
ATTAGGCTCCGTAACCCGCAAG
+
!!!!!!!!!!!!!!!!!!!!!!
@SEQ_ID3
GCCACCCCTTGGTCTTTATGCA
+
IIIIIIIIIIIIIIIIIIIIII
EOF

cat > "$dir_in/b_2.fq" <<'EOF'
@SEQ_ID1
AGCCAGCAGCAATGTGGTAAAG
+
IIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
CTTGCGGGTTACGGAGCCTAAT
+
!!!!!!!!!!!!!!!!!!!!!!
@SEQ_ID3
TGCATAAAGACCAAGGGGTGGC
+
IIIIIIIIIIIIIIIIIIIIII
EOF

echo "==============================================================================="
echo "> Run salmon index"

salmon index \
  --transcripts "$dir_in/transcriptome.fasta" \
  --index "$dir_in/index" \
  --kmerLen 11

echo "==============================================================================="
echo "> Run salmon quant for single-end reads"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/index" \
  --unmated_reads "$dir_in/a_1.fq" \
  --output "quant_se_results" \
  --quant_results "quant_se.sf" \
  --min_assigned_frags 1

echo ">> Checking output"
[ ! -d "quant_se_results" ] && echo "Output directory quant_se_results does not exist" && exit 1
[ ! -f "quant_se.sf" ] && echo "Output file quant_se.sf does not exist!" && exit 1
[ ! -s "quant_se.sf" ] && echo "Output file quant_se.sf is empty!" && exit 1
grep -q "Name	Length	EffectiveLength	TPM	NumReads" "quant_se.sf" || (echo "Output file quant_se.sf does not have the right format!" && exit 1)
[ $(grep "contig1" "quant_se.sf" | cut -f 5) != '2.000' ] && echo "Number of reads mapping to contig1 does not match the expected value!" && exit 1
[ $(grep "contig2" "quant_se.sf" | cut -f 5) != '0.000' ] && echo "Number of reads mapping to contig2 does not match the expected value!" && exit 1
[ $(grep '"percent_mapped":' quant_se_results/aux_info/meta_info.json | cut -d':' -f 2) != '66.66666666666666,' ] && echo "Mapping rate does not match the expected value!" && exit 1

echo "==============================================================================="
echo "> Run salmon quant for paired-end reads"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/index" \
  --mates1 "$dir_in/a_1.fq" \
  --mates2 "$dir_in/a_2.fq" \
  --output "quant_pe_results" \
  --quant_results "quant_pe.sf" \
  --min_assigned_frags 1

echo ">> Checking output"
[ ! -d "quant_pe_results" ] && echo "Output directory quant_pe_results does not exist" && exit 
[ ! -f "quant_pe.sf" ] && echo "Output file quant_pe.sf does not exist!" && exit 1
[ ! -s "quant_pe.sf" ] && echo "Output file quant_pe.sf is empty!" && exit 1
grep -q "Name	Length	EffectiveLength	TPM	NumReads" "quant_pe.sf" || (echo "Output file quant_pe.sf does not have the right format!" && exit 1)
[ $(grep "contig1" "quant_pe.sf" | cut -f 5) != '2.000' ] && echo "Number of reads mapping to contig1 does not match the expected value!" && exit 1
[ $(grep "contig2" "quant_pe.sf" | cut -f 5) != '0.000' ] && echo "Number of reads mapping to contig2 does not match the expected value!" && exit 1
[ $(grep '"percent_mapped":' quant_pe_results/aux_info/meta_info.json | cut -d':' -f 2) != '66.66666666666666,' ] && echo "Mapping rate does not match the expected value!" && exit 1

echo "==============================================================================="
echo "> Run salmon quant for paired-end reads with technical replicates"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/index" \
  --mates1 "$dir_in/a_1.fq;$dir_in/b_1.fq" \
  --mates2 "$dir_in/a_2.fq;$dir_in/b_2.fq" \
  --output "quant_pe_rep_results" \
  --quant_results "quant_pe_rep.sf" \
  --min_assigned_frags 1

echo ">> Checking output"
[ ! -d "quant_pe_rep_results" ] && echo "Output directory quant_pe_rep_results does not exist" && exit 
[ ! -f "quant_pe_rep.sf" ] && echo "Output file quant_pe_rep.sf does not exist!" && exit 1
[ ! -s "quant_pe_rep.sf" ] && echo "Output file quant_pe_rep.sf is empty!" && exit 1
grep -q "Name	Length	EffectiveLength	TPM	NumReads" "quant_pe_rep.sf" || (echo "Output file quant_pe_rep.sf does not have the right format!" && exit 1)
[ $(grep "contig1" "quant_pe_rep.sf" | cut -f 5) != '2.000' ] && echo "Number of reads mapping to contig1 does not match the expected value!" && exit 1
[ $(grep "contig2" "quant_pe_rep.sf" | cut -f 5) != '2.000' ] && echo "Number of reads mapping to contig2 does not match the expected value!" && exit 1
[ $(grep '"percent_mapped":' quant_pe_rep_results/aux_info/meta_info.json | cut -d':' -f 2) != '66.66666666666666,' ] && echo "Mapping rate does not match the expected value!" && exit 1


# TODO: check counts and mapping rates
# contig1 should have 2 reads, contig2 should have 2 reads
# mapping rate should be 66.6%

echo "==============================================================================="
echo "> Test successful"