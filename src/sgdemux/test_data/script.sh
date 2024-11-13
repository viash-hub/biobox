#!/bin/bash
set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

OUT=src/sgdemux/test_data/


TAR_LOC="$OUT/unfiltered_fastq.tar"
if [ ! -f "$TAR_LOC" ]; then
    wget https://singular-public-repo.s3.us-west-1.amazonaws.com/example_raw_files/unfiltered_fastq.tar.gz -O "$TAR_LOC"
fi

tar -xvf "$TAR_LOC" -C "$OUT"

# NOTE: sgdemux requires block compressed gzip files!
function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing $(basename $input)"
    seqkit head -n 10000 "$input" | bgzip --threads 12 > "$output"
  fi
}
tar_contents=(
    Undetermined_S0_L001_I1_001.fastq.gz
    Undetermined_S0_L001_I2_001.fastq.gz
    Undetermined_S0_L001_R1_001.fastq.gz
    Undetermined_S0_L001_R2_001.fastq.gz
    Undetermined_S0_L002_I1_001.fastq.gz
    Undetermined_S0_L002_I2_001.fastq.gz
    Undetermined_S0_L002_R1_001.fastq.gz
    Undetermined_S0_L002_R2_001.fastq.gz
    Undetermined_S0_L003_I1_001.fastq.gz
    Undetermined_S0_L003_I2_001.fastq.gz
    Undetermined_S0_L003_R1_001.fastq.gz
    Undetermined_S0_L003_R2_001.fastq.gz
    Undetermined_S0_L004_I1_001.fastq.gz
    Undetermined_S0_L004_I2_001.fastq.gz
    Undetermined_S0_L004_R1_001.fastq.gz
    Undetermined_S0_L004_R2_001.fastq.gz
)

mkdir -p "$OUT/fastq"
for fastq in ${tar_contents[@]}; do
    seqkit_head "$OUT/unfiltered_fastq/$fastq" "$OUT/fastq/$fastq"
done
cp "$OUT/unfiltered_fastq/samplesheet.csv" "$OUT/samplesheet.csv"