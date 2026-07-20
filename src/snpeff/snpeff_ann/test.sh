#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

###########################################################################

# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
log "Created $TMPDIR"
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT
DATA_DIR="$TMPDIR/data"
mkdir "$DATA_DIR"
TEST_GENOME="$DATA_DIR/test"
mkdir "$TEST_GENOME"
log "Downloading test data"
curl -o - https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz > "$TEST_GENOME/sequences.fa.gz"
curl -o - https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz | zcat | sed -n '/^[#1]\s/p' | gzip > "$TEST_GENOME/genes.gtf.gz"
snpEff build -dataDir "$DATA_DIR" -noCheckCds -noCheckProtein -gtf22 -noLog -configOption "test.genome=GRCh38Chr1" -v test

# Test 1: Run SnpEff with only required parameters

log "> Run Test 1: required parameters"
mkdir test1
pushd test1 > /dev/null # cd test1 (stack)

"$meta_executable" \
  --genome_version test \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --input "$meta_resources_dir/test_data/cancer.vcf" \
  --output out.vcf

check_file_exists "out.vcf" "annotated VCF output"
check_file_not_empty "out.vcf" "annotated VCF output"

popd > /dev/null # Remove directory from stack (LIFO)

log "Test 1 succeeded."

###########################################################################

# Test 2: Run SnpEff with a different input + options

log "> Run Test 2: different input + options"
mkdir test2
pushd test2 > /dev/null

"$meta_executable" \
  --genome_version test \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --input "$meta_resources_dir/test_data/test.vcf" \
  --interval "$meta_resources_dir/test_data/my_annotations.bed" \
  --no_stats \
  --output output.vcf

check_file_exists "output.vcf" "annotated VCF output"
check_file_not_empty "output.vcf" "annotated VCF output"

check_file_not_exists "snpEff_summary.html" "HTML summary (disabled by --no_stats)"
check_file_not_exists "snpEff_summary.genes.txt" "genes stats file (disabled by --no_stats)"

popd > /dev/null

log "Test 2 succeeded."

###########################################################################

# Test 3: Move the output files to other locations

log "> Run Test 3: move output files"
mkdir test3
pushd test3 > /dev/null

"$meta_executable" \
  --genome_version test \
  --input "$meta_resources_dir/test_data/test.vcf" \
  --output output.vcf \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --summary temp/summary.html \
  --genes temp/genes.txt

check_file_exists "output.vcf" "annotated VCF output"
check_file_not_empty "output.vcf" "annotated VCF output"

check_file_exists "temp/summary.html" "HTML summary (relocated via --summary)"
check_file_not_empty "temp/summary.html" "HTML summary (relocated via --summary)"

check_file_exists "temp/genes.txt" "genes stats file (relocated via --genes)"
check_file_not_empty "temp/genes.txt" "genes stats file (relocated via --genes)"

popd > /dev/null

log "Test 3 succeeded."

###########################################################################

# Test 4: --stats names the intermediate file, --summary sets the final path

log "> Run Test 4: --stats intermediate name + --summary final path"
mkdir test4
pushd test4 > /dev/null

"$meta_executable" \
  --genome_version test \
  --input "$meta_resources_dir/test_data/test.vcf" \
  --output output.vcf \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --stats custom_intermediate.html \
  --summary final_summary.html

check_file_exists "output.vcf" "annotated VCF output"
check_file_not_empty "output.vcf" "annotated VCF output"

# The HTML summary should end up at the --summary path, not the --stats path.
check_file_exists "final_summary.html" "HTML summary (final --summary path)"
check_file_not_empty "final_summary.html" "HTML summary (final --summary path)"
check_file_not_exists "custom_intermediate.html" "intermediate --stats file (should have been moved)"

popd > /dev/null

log "Test 4 succeeded."

###########################################################################

# Test 5: --no_hgvs, --only_tr and --fastaprot/--fastaprot_no_ref

log "> Run Test 5: --no_hgvs, --only_tr, --fastaprot_no_ref"
mkdir test5
pushd test5 > /dev/null

# ENST00000641515 (OR4F5) is a real, protein-coding chr1 transcript whose CDS
# overlaps the variant in test_data/cancer.vcf (position 69091)
echo "ENST00000641515" > only_tr.txt

"$meta_executable" \
  --genome_version test \
  --input "$meta_resources_dir/test_data/cancer.vcf" \
  --output output.vcf \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --no_hgvs \
  --only_tr only_tr.txt \
  --fastaprot proteins.fa \
  --fastaprot_no_ref

check_file_exists "output.vcf" "annotated VCF output"
check_file_not_empty "output.vcf" "annotated VCF output"

# --only_tr should restrict annotation to just the listed transcript
check_file_contains "output.vcf" "ENST00000641515" "annotated VCF output (restricted via --only_tr)"

check_file_exists "proteins.fa" "protein FASTA"
check_file_not_empty "proteins.fa" "protein FASTA"

popd > /dev/null

log "Test 5 succeeded."

###########################################################################

print_test_summary "snpeff_ann"
