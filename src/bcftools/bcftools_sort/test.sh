#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

test_data="$meta_resources_dir/test_data"

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_identical_content() {
  diff -a "$2" "$1" \
    || (echo "Files are not identical!" && exit 1)
}
#############################################

# Create directories for tests
echo "Creating Test Data..."
TMPDIR=$(mktemp -d "$meta_temp_dir/XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Create test data
cat <<EOF > "$TMPDIR/example.vcf"
##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=58617616>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
##ALT=<ID=CNV,Description="Copy number variable region">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
19	112	.	A	G	10	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3
19	111	.	A	C	9.6	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3
20	1235237	.	T	.	.	.	.	GT	0/0	0|0	./.
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.,.
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4:.,.
20	1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:.:56,60	0|0:48:4:51,51	0/0:61:2:.,.
20	1234567	microsat1	G	GA,GAC	50	PASS	NS=3;DP=9;AA=G;AN=6;AC=3,1	GT:GQ:DP	0/1:.:4	0/2:17:2	1/1:40:3
EOF

# Create expected output
cat <<EOF > "$TMPDIR/expected_output.vcf"
##fileformat=VCFv4.0
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=58617616>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
##ALT=<ID=CNV,Description="Copy number variable region">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
19	111	.	A	C	9.6	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3
19	112	.	A	G	10	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.,.
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4:.,.
20	1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:.:56,60	0|0:48:4:51,51	0/0:61:2:.,.
20	1234567	microsat1	G	GA,GAC	50	PASS	NS=3;DP=9;AA=G;AN=6;AC=3,1	GT:GQ:DP	0/1:.:4	0/2:17:2	1/1:40:3
20	1235237	.	T	.	.	.	.	GT	0/0	0|0	./.
EOF

cat <<EOF > "$TMPDIR/expected_bcf.vcf"
##fileformat=VCFv4.0
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=58617616>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
##ALT=<ID=CNV,Description="Copy number variable region">
##bcftools_viewVersion=1.16+htslib-1.16
##bcftools_viewCommand=view -O b -o example.bcf example.vcf.gz; Date=Mon Aug 26 13:00:22 2024
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
19	111	.	A	C	9.6	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3
19	112	.	A	G	10	.	.	GT:HQ	0|0:10,10	0|0:10,10	0/1:3,3
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.,.
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4:.,.
20	1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:.:56,60	0|0:48:4:51,51	0/0:61:2:.,.
20	1234567	microsat1	G	GA,GAC	50	PASS	NS=3;DP=9;AA=G;AN=6;AC=3,1	GT:GQ:DP	0/1:.:4	0/2:17:2	1/1:40:3
20	1235237	.	T	.	.	.	.	GT	0/0	0|0	./.
EOF


# Test 1: Default Use
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bcftools_sort on VCF file"
"$meta_executable" \
  --input "../example.vcf" \
  --output "output.vcf" \
  --output_type "v" \
  &> /dev/null

# checks
assert_file_exists "output.vcf"
assert_file_not_empty "output.vcf"
assert_identical_content "output.vcf" "../expected_output.vcf"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: BCF file input
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bcftools_sort on BCF file"
"$meta_executable" \
  --input "${test_data}/example.bcf" \
  --output "output.vcf" \
  --output_type "v" \
  &> /dev/null

# checks
assert_file_exists "output.vcf"
assert_file_not_empty "output.vcf"
assert_identical_content "output.vcf" "../expected_bcf.vcf"
echo "- test2 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
