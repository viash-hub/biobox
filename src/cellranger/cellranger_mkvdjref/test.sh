#!/bin/bash

## VIASH START
## VIASH END

# source the helpers & setup test env
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Test for $meta_name starting ..."

# count the number of records in a FASTA file
count_fasta_records() {
  grep -c "^>" "$1" || true
}

# create test data: a small contig with four annotated V(D)J genes
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"
create_test_fasta "$test_data_dir/genome.fa" 1 6000
input_gtf="$test_data_dir/genes.gtf"

cat > "$input_gtf" << 'EOF'
seq1	test	gene	1000	1300	.	+	.	gene_id "ENSG01"; gene_name "TRAV1-1"; gene_biotype "TR_V_gene";
seq1	test	CDS	1000	1300	.	+	0	gene_id "ENSG01"; transcript_id "ENST01"; gene_name "TRAV1-1"; gene_biotype "TR_V_gene"; transcript_biotype "TR_V_gene";
seq1	test	gene	2000	2060	.	+	.	gene_id "ENSG02"; gene_name "TRAJ1"; gene_biotype "TR_J_gene";
seq1	test	CDS	2000	2060	.	+	0	gene_id "ENSG02"; transcript_id "ENST02"; gene_name "TRAJ1"; gene_biotype "TR_J_gene"; transcript_biotype "TR_J_gene";
seq1	test	gene	3000	3500	.	+	.	gene_id "ENSG03"; gene_name "TRAC"; gene_biotype "TR_C_gene";
seq1	test	CDS	3000	3500	.	+	0	gene_id "ENSG03"; transcript_id "ENST03"; gene_name "TRAC"; gene_biotype "TR_C_gene"; transcript_biotype "TR_C_gene";
seq1	test	gene	4000	4300	.	-	.	gene_id "ENSG04"; gene_name "TRBV1"; gene_biotype "TR_V_gene";
seq1	test	CDS	4000	4300	.	-	0	gene_id "ENSG04"; transcript_id "ENST04"; gene_name "TRBV1"; gene_biotype "TR_V_gene"; transcript_biotype "TR_V_gene";
EOF

check_file_exists "$input_gtf" "test input GTF"

# --- TEST_1: build a reference from a genome FASTA and a gene GTF ---
log "Starting TEST_1: build a V(D)J reference from --fasta and --genes"

output_dir="$meta_temp_dir/test1/reference"
"$meta_executable" \
  --genome test_vdj_ref \
  --fasta "$test_data_dir/genome.fa" \
  --genes "$input_gtf" \
  --output "$output_dir"

check_dir_exists "$output_dir" "V(D)J reference"
check_file_exists "$output_dir/fasta/regions.fa" "V(D)J regions FASTA"
check_file_exists "$output_dir/reference.json" "V(D)J reference metadata"
check_file_contains "$output_dir/reference.json" '"genomes": "test_vdj_ref"' "reference metadata"

num_records=$(count_fasta_records "$output_dir/fasta/regions.fa")
[[ "$num_records" -eq 4 ]] || { log_error "Expected 4 V(D)J regions, found $num_records"; exit 1; }
check_file_contains "$output_dir/fasta/regions.fa" "TRAV1-1|L-REGION+V-REGION|TR|TRA" "V(D)J regions FASTA"
check_file_contains "$output_dir/fasta/regions.fa" "TRAJ1|J-REGION|TR|TRA" "V(D)J regions FASTA"
check_file_contains "$output_dir/fasta/regions.fa" "TRAC|C-REGION|TR|TRA" "V(D)J regions FASTA"
check_file_contains "$output_dir/fasta/regions.fa" "TRBV1|L-REGION+V-REGION|TR|TRB" "V(D)J regions FASTA"

log "✓ TEST_1 completed successfully"

# --- TEST_2: build a reference from a V(D)J segment FASTA ---
log "Starting TEST_2: build a V(D)J reference from --seqs"

# the regions FASTA written by TEST_1 already follows the mkvdjref spec
segments_fa="$test_data_dir/segments.fa"
cp "$output_dir/fasta/regions.fa" "$segments_fa"

output_dir="$meta_temp_dir/test2/reference"
"$meta_executable" \
  --genome test_vdj_seqs_ref \
  --seqs "$segments_fa" \
  --output "$output_dir"

check_file_exists "$output_dir/fasta/regions.fa" "V(D)J regions FASTA"
check_file_contains "$output_dir/reference.json" '"input_gtf_files": null' "reference metadata"

num_records=$(count_fasta_records "$output_dir/fasta/regions.fa")
[[ "$num_records" -eq 4 ]] || { log_error "Expected 4 V(D)J regions, found $num_records"; exit 1; }

log "✓ TEST_2 completed successfully"

# --- TEST_3: drop transcripts and record a reference version ---
log "Starting TEST_3: --rm_transcripts and --ref_version"

rm_transcripts="$test_data_dir/remove_transcripts.txt"
echo "ENST04" > "$rm_transcripts"

output_dir="$meta_temp_dir/test3/reference"
"$meta_executable" \
  --genome test_vdj_ref \
  --fasta "$test_data_dir/genome.fa" \
  --genes "$input_gtf" \
  --rm_transcripts "$rm_transcripts" \
  --ref_version "1.2.3" \
  --output "$output_dir"

check_file_contains "$output_dir/reference.json" '"version": "1.2.3"' "reference metadata"
check_file_not_contains "$output_dir/fasta/regions.fa" "TRBV1" "V(D)J regions FASTA"

num_records=$(count_fasta_records "$output_dir/fasta/regions.fa")
[[ "$num_records" -eq 3 ]] || { log_error "Expected 3 V(D)J regions, found $num_records"; exit 1; }

log "✓ TEST_3 completed successfully"

# --- TEST_4: gzipped inputs are decompressed ---
log "Starting TEST_4: gzipped --fasta and --genes"

gzip -c "$test_data_dir/genome.fa" > "$test_data_dir/genome.fa.gz"
gzip -c "$input_gtf" > "$test_data_dir/genes.gtf.gz"

output_dir="$meta_temp_dir/test4/reference"
"$meta_executable" \
  --genome test_vdj_ref \
  --fasta "$test_data_dir/genome.fa.gz" \
  --genes "$test_data_dir/genes.gtf.gz" \
  --output "$output_dir"

num_records=$(count_fasta_records "$output_dir/fasta/regions.fa")
[[ "$num_records" -eq 4 ]] || { log_error "Expected 4 V(D)J regions, found $num_records"; exit 1; }

log "✓ TEST_4 completed successfully"

# --- TEST_5: --seqs and --fasta are mutually exclusive ---
log "Starting TEST_5: --seqs combined with --fasta is rejected"

output_dir="$meta_temp_dir/test5/reference"
if "$meta_executable" \
  --genome test_vdj_ref \
  --fasta "$test_data_dir/genome.fa" \
  --seqs "$segments_fa" \
  --output "$output_dir" > /dev/null 2>&1; then
  log_error "✗ Expected a non-zero exit code when combining --seqs with --fasta"
  exit 1
fi
log "✓ Combining --seqs with --fasta failed as expected"

log "✓ TEST_5 completed successfully"

print_test_summary "$meta_name tests"
