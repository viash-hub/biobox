#!/bin/bash

## VIASH START
## VIASH END

# source the helpers & setup test env
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Test for $meta_name starting ..."

# create test data: three genes, each with a different gene_biotype
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"
input_gtf="$test_data_dir/genes.gtf"

cat > "$input_gtf" << 'EOF'
chr1	test	gene	1000	1999	.	+	.	gene_id "gene1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1	test	transcript	1000	1999	.	+	.	gene_id "gene1"; transcript_id "transcript1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1	test	exon	1000	1999	.	+	.	gene_id "gene1"; transcript_id "transcript1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1	test	gene	3000	3999	.	+	.	gene_id "gene2"; gene_name "GENE2"; gene_biotype "lincRNA";
chr1	test	transcript	3000	3999	.	+	.	gene_id "gene2"; transcript_id "transcript2"; gene_name "GENE2"; gene_biotype "lincRNA";
chr1	test	exon	3000	3999	.	+	.	gene_id "gene2"; transcript_id "transcript2"; gene_name "GENE2"; gene_biotype "lincRNA";
chr1	test	gene	5000	5999	.	-	.	gene_id "gene3"; gene_name "GENE3"; gene_biotype "pseudogene";
chr1	test	transcript	5000	5999	.	-	.	gene_id "gene3"; transcript_id "transcript3"; gene_name "GENE3"; gene_biotype "pseudogene";
chr1	test	exon	5000	5999	.	-	.	gene_id "gene3"; transcript_id "transcript3"; gene_name "GENE3"; gene_biotype "pseudogene";
EOF

check_file_exists "$input_gtf" "test input GTF"

# --- TEST_1: filter on a single attribute ---
log "Starting TEST_1: filter on a single key-value pair"

output_gtf="$meta_temp_dir/test1/genes_filtered.gtf"
"$meta_executable" \
  --input_gtf "$input_gtf" \
  --output_gtf "$output_gtf" \
  --attribute "gene_biotype:protein_coding"

check_file_exists "$output_gtf" "filtered GTF"
check_file_line_count "$output_gtf" 3 "filtered GTF"
check_file_contains "$output_gtf" 'gene_id "gene1"' "filtered GTF"
check_file_not_contains "$output_gtf" 'gene_id "gene2"' "filtered GTF"
check_file_not_contains "$output_gtf" 'gene_id "gene3"' "filtered GTF"

log "✓ TEST_1 completed successfully"

# --- TEST_2: filter on multiple values for the same attribute key ---
log "Starting TEST_2: filter on multiple values for the same key"

output_gtf="$meta_temp_dir/test2/genes_filtered.gtf"
"$meta_executable" \
  --input_gtf "$input_gtf" \
  --output_gtf "$output_gtf" \
  --attribute "gene_biotype:protein_coding" \
  --attribute "gene_biotype:lincRNA"

check_file_exists "$output_gtf" "filtered GTF"
check_file_line_count "$output_gtf" 6 "filtered GTF"
check_file_contains "$output_gtf" 'gene_id "gene1"' "filtered GTF"
check_file_contains "$output_gtf" 'gene_id "gene2"' "filtered GTF"
check_file_not_contains "$output_gtf" 'gene_id "gene3"' "filtered GTF"

log "✓ TEST_2 completed successfully"

# --- TEST_3: no attribute filtering ---
log "Starting TEST_3: no --attribute, all records are retained"

output_gtf="$meta_temp_dir/test3/genes_filtered.gtf"
"$meta_executable" \
  --input_gtf "$input_gtf" \
  --output_gtf "$output_gtf"

check_file_exists "$output_gtf" "filtered GTF"
check_file_line_count "$output_gtf" 9 "filtered GTF"
check_file_contains "$output_gtf" 'gene_id "gene3"' "filtered GTF"

log "✓ TEST_3 completed successfully"

print_test_summary "$meta_name tests"
