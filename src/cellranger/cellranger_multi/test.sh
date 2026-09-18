#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting tests for $meta_name"

log "Copying test data from the Cell Ranger installation directory"
test_data="$meta_temp_dir/test_data"
mkdir -p "$test_data"
cp -r "/opt/cellranger-10.0.0/external/cellranger_tiny_fastq" "$test_data"
cp -r "/opt/cellranger-10.0.0/external/cellranger_tiny_ref" "$test_data"
fastq_dir="$test_data/cellranger_tiny_fastq"
reference_dir="$test_data/cellranger_tiny_ref"

log "Creating a gzipped tarball of the reference"
reference_tar="$test_data/cellranger_tiny_ref.tar.gz"
tar -czf "$reference_tar" -C "$test_data" "cellranger_tiny_ref"

##############################################################
log "Starting TEST 1: staging into an existing multi config CSV"
##############################################################
# The paths in the CSV are placeholders: they are replaced by the staged locations of
# --gex_reference and --fastqs.
cat > "$meta_temp_dir/test1_config.csv" <<'EOF'
[gene-expression]
reference,/placeholder/transcriptome
create-bam,false

[libraries]
fastq_id,fastqs,feature_types,lanes
tinygex,/placeholder/fastqs,Gene Expression,1
EOF

output_dir="$meta_temp_dir/test1"
"$meta_executable" \
  --csv "$meta_temp_dir/test1_config.csv" \
  --fastqs "$fastq_dir" \
  --gex_reference "$reference_dir" \
  --id "tinygex" \
  --output "$output_dir"

check_file_exists "$output_dir/config.csv" "Multi config CSV used for the run"
check_file_not_contains "$output_dir/config.csv" "/placeholder" "Multi config CSV used for the run"
check_file_exists "$output_dir/raw_feature_bc_matrix.h5" "Raw feature-barcode matrix"
check_dir_exists "$output_dir/per_sample_outs/tinygex" "Per sample output directory"
check_file_exists "$output_dir/per_sample_outs/tinygex/sample_filtered_feature_bc_matrix.h5" \
  "Sample filtered feature-barcode matrix"
check_file_exists "$output_dir/per_sample_outs/tinygex/metrics_summary.csv" "Sample metrics summary"
check_file_contains "$output_dir/cellranger_multi.log" "Pipestance completed successfully" \
  "Cell Ranger log"
log "TEST 1 completed successfully"

##############################################################
log "Starting TEST 2: appending an entry, a tarballed reference and individual FASTQ files"
##############################################################
# The [gene-expression] section has no reference entry, so it is appended to the section.
cat > "$meta_temp_dir/test2_config.csv" <<'EOF'
[gene-expression]
create-bam,false

[libraries]
fastq_id,fastqs,feature_types,lanes
tinygex,/placeholder/fastqs,Gene Expression,1
EOF

output_dir="$meta_temp_dir/test2"
"$meta_executable" \
  --csv "$meta_temp_dir/test2_config.csv" \
  --fastqs "$fastq_dir/tinygex_S1_L001_R1_001.fastq.gz" \
  --fastqs "$fastq_dir/tinygex_S1_L001_R2_001.fastq.gz" \
  --gex_reference "$reference_tar" \
  --description "A tiny gene expression run" \
  --output "$output_dir"

check_file_contains "$output_dir/config.csv" "^reference,/" "Multi config CSV used for the run"
check_file_not_contains "$output_dir/config.csv" "/placeholder" "Multi config CSV used for the run"
check_file_exists "$output_dir/per_sample_outs/run/sample_filtered_feature_bc_matrix.h5" \
  "Sample filtered feature-barcode matrix"
check_file_contains "$output_dir/per_sample_outs/run/web_summary.html" "A tiny gene expression run" \
  "Web summary"
log "TEST 2 completed successfully"

##############################################################
log "Starting TEST 3: dry run, appending a section that the CSV does not have"
##############################################################
# Cell Ranger rejects a CSV that declares VDJ libraries without a [vdj] section, so the run
# only gets as far as generating the pipeline invocation if --vdj_reference was appended as
# a new section. A dry run does not read the reference itself, so the gene expression
# reference stands in for a V(D)J one here.
cat > "$meta_temp_dir/test3_config.csv" <<'EOF'
[libraries]
fastq_id,fastqs,feature_types,lanes
tinyvdj,/placeholder/fastqs,VDJ-T,1
EOF

output_dir="$meta_temp_dir/test3"
"$meta_executable" \
  --csv "$meta_temp_dir/test3_config.csv" \
  --fastqs "$fastq_dir" \
  --vdj_reference "$reference_dir" \
  --dry \
  --output "$output_dir"

check_dir_exists "$output_dir" "Output directory"
check_file_contains "$output_dir/config.csv" "^\[vdj\]" "Multi config CSV generated for the run"
check_file_contains "$output_dir/config.csv" "^reference,/" "Multi config CSV generated for the run"
check_file_not_contains "$output_dir/config.csv" "/placeholder" "Multi config CSV generated for the run"
log "TEST 3 completed successfully"

##############################################################
log "Starting TEST 4: multi config CSV without a fastqs column"
##############################################################
cat > "$meta_temp_dir/test4_config.csv" <<'EOF'
[gene-expression]
reference,/placeholder/transcriptome
create-bam,false

[libraries]
fastq_id,feature_types,lanes
tinygex,Gene Expression,1
EOF

output_dir="$meta_temp_dir/test4"
if "$meta_executable" \
  --csv "$meta_temp_dir/test4_config.csv" \
  --fastqs "$fastq_dir" \
  --gex_reference "$reference_dir" \
  --output "$output_dir" > "$meta_temp_dir/test4.log" 2>&1; then
  log_error "Expected the component to fail on a [libraries] section without a fastqs column"
  exit 1
fi
check_file_contains "$meta_temp_dir/test4.log" "has no fastqs column" "Error output"
log "TEST 4 completed successfully"

##############################################################
log "Starting TEST 4b: a run that Cell Ranger itself rejects keeps the log"
##############################################################
# create-bam is a required entry, so Cell Ranger exits before doing any work. The log is
# written while the pipeline runs, so it has to survive the failure.
cat > "$meta_temp_dir/test4b_config.csv" <<'EOF'
[libraries]
fastq_id,fastqs,feature_types,lanes
tinygex,/placeholder/fastqs,Gene Expression,1
EOF

output_dir="$meta_temp_dir/test4b"
if "$meta_executable" \
  --csv "$meta_temp_dir/test4b_config.csv" \
  --fastqs "$fastq_dir" \
  --gex_reference "$reference_dir" \
  --output "$output_dir" > "$meta_temp_dir/test4b.log" 2>&1; then
  log_error "Expected the component to fail on a config without create-bam"
  exit 1
fi
check_file_contains "$output_dir/cellranger_multi.log" "create-bam is a required parameter" \
  "Cell Ranger log of the failed run"
log "TEST 4b completed successfully"

##############################################################
log "Starting TEST 5: the same library sequenced on two flowcells"
##############################################################
# A library sequenced on more than one flowcell yields identically named FASTQ files in
# different directories, which only the fastqs column tells apart. The directories that hold
# them are both called "a", so they can only be told apart by more than the last component of
# their path. The second flowcell holds the lane 2 reads under the file names of lane 1,
# uncompressed, which also covers that Cell Ranger reads uncompressed FASTQ files.
mkdir -p "$test_data/flowcell1/a" "$test_data/flowcell2/a"
cp "$fastq_dir/tinygex_S1_L001_R1_001.fastq.gz" "$test_data/flowcell1/a/"
cp "$fastq_dir/tinygex_S1_L001_R2_001.fastq.gz" "$test_data/flowcell1/a/"
gunzip -c "$fastq_dir/tinygex_S1_L002_R1_001.fastq.gz" > "$test_data/flowcell2/a/tinygex_S1_L001_R1_001.fastq"
gunzip -c "$fastq_dir/tinygex_S1_L002_R2_001.fastq.gz" > "$test_data/flowcell2/a/tinygex_S1_L001_R2_001.fastq"

cat > "$meta_temp_dir/test5_config.csv" <<'EOF'
[gene-expression]
create-bam,false

[libraries]
fastq_id,fastqs,feature_types
tinygex,/placeholder/flowcell1/a,Gene Expression
tinygex,/placeholder/flowcell2/a,Gene Expression
EOF

output_dir="$meta_temp_dir/test5"
"$meta_executable" \
  --csv "$meta_temp_dir/test5_config.csv" \
  --fastqs "$test_data/flowcell1/a" \
  --fastqs "$test_data/flowcell2/a" \
  --gex_reference "$reference_dir" \
  --output "$output_dir"

# Each row keeps pointing at the reads it came from, so the two rows end up on two different
# staged directories rather than on one.
check_file_contains "$output_dir/config.csv" "/a,Gene Expression" "Multi config CSV used for the run"
check_file_contains "$output_dir/config.csv" "/a_2,Gene Expression" "Multi config CSV used for the run"
# Cell Ranger reports the reads of the second directory separately, as "tinygex (2)".
check_file_contains "$output_dir/per_sample_outs/run/metrics_summary.csv" "tinygex (2)" \
  "Sample metrics summary"
log "TEST 5 completed successfully"

##############################################################
log "Starting TEST 6: a fastqs column that does not identify one input directory"
##############################################################
# The value matches none of the input directories, and then matches both equally well.
run_expecting_failure() {
  local name="$1" fastqs_value="$2"
  cat > "$meta_temp_dir/${name}_config.csv" <<EOF
[gene-expression]
create-bam,false

[libraries]
fastq_id,fastqs,feature_types
tinygex,${fastqs_value},Gene Expression
EOF
  if "$meta_executable" \
    --csv "$meta_temp_dir/${name}_config.csv" \
    --fastqs "$test_data/flowcell1/a" \
    --fastqs "$test_data/flowcell2/a" \
    --gex_reference "$reference_dir" \
    --output "$meta_temp_dir/$name" > "$meta_temp_dir/$name.log" 2>&1; then
    log_error "Expected the component to fail on a fastqs column of '${fastqs_value}'"
    exit 1
  fi
  check_file_contains "$meta_temp_dir/$name.log" \
    "does not match exactly one of the directories the FASTQ files came from" "Error output"
}

run_expecting_failure "test6a" "/placeholder/somewhere_else"
run_expecting_failure "test6b" "a"
log "TEST 6 completed successfully"

##############################################################
log "Starting TEST 7: lz4 compressed FASTQ files"
##############################################################
# Cell Ranger accepts .fastq, .fastq.gz and .fastq.lz4. A dry run is enough here: staging
# happens before Cell Ranger is called, so the run only gets this far if the lz4 files were
# picked up at all.
mkdir -p "$test_data/lz4"
lz4 -q "$test_data/flowcell2/a/tinygex_S1_L001_R1_001.fastq" "$test_data/lz4/tinygex_S1_L001_R1_001.fastq.lz4"
lz4 -q "$test_data/flowcell2/a/tinygex_S1_L001_R2_001.fastq" "$test_data/lz4/tinygex_S1_L001_R2_001.fastq.lz4"

output_dir="$meta_temp_dir/test7"
"$meta_executable" \
  --csv "$meta_temp_dir/test1_config.csv" \
  --fastqs "$test_data/lz4" \
  --gex_reference "$reference_dir" \
  --dry \
  --output "$output_dir"

check_file_contains "$output_dir/config.csv" "/lz4,Gene Expression" "Multi config CSV generated for the run"
log "TEST 7 completed successfully"

print_test_summary "All tests"
