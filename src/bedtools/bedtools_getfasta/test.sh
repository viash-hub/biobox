#!/usr/bin/env bash
set -eo pipefail

TMPDIR=$(mktemp -d)
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Create dummy test fasta file
cat > "$TMPDIR/test.fa" <<EOF
>chr1
AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG
EOF

TAB="$(printf '\t')"

# Create dummy bed file
cat > "$TMPDIR/test.bed" <<EOF
chr1${TAB}5${TAB}10${TAB}myseq
EOF

# Create expected bed file
cat > "$TMPDIR/expected.fasta" <<EOF
>chr1:5-10
AAACC
EOF

"$meta_executable" \
  --input_bed "$TMPDIR/test.bed" \
  --input_fasta "$TMPDIR/test.fa" \
  --output "$TMPDIR/output.fasta"

cmp --silent "$TMPDIR/output.fasta" "$TMPDIR/expected.fasta" || { echo "files are different:"; exit 1; }


# Create expected bed file for --name
cat > "$TMPDIR/expected_with_name.fasta" <<EOF
>myseq::chr1:5-10
AAACC
EOF

"$meta_executable" \
  --input_bed "$TMPDIR/test.bed" \
  --input_fasta "$TMPDIR/test.fa" \
  --name \
  --output "$TMPDIR/output_with_name.fasta"


cmp --silent "$TMPDIR/output_with_name.fasta" "$TMPDIR/expected_with_name.fasta" || { echo "Files when using --name are different."; exit 1; }

# Create expected bed file for --name_only
cat > "$TMPDIR/expected_with_name_only.fasta" <<EOF
>myseq
AAACC
EOF

"$meta_executable" \
  --input_bed "$TMPDIR/test.bed" \
  --input_fasta "$TMPDIR/test.fa" \
  --name_only \
  --output "$TMPDIR/output_with_name_only.fasta"

cmp --silent "$TMPDIR/output_with_name_only.fasta" "$TMPDIR/expected_with_name_only.fasta" || { echo "Files when using --name_only are different."; exit 1; }


# Create expected tab-delimited file for --tab
cat > "$TMPDIR/expected_tab.out" <<EOF
myseq${TAB}AAACC
EOF

"$meta_executable" \
  --input_bed "$TMPDIR/test.bed" \
  --input_fasta "$TMPDIR/test.fa" \
  --name_only \
  --tab \
  --output "$TMPDIR/tab.out"

cmp --silent "$TMPDIR/expected_tab.out" "$TMPDIR/tab.out" || { echo "Files when using --tab are different."; exit 1; }


# Create expected tab-delimited file for --bed_out
cat > "$TMPDIR/expected.bed" <<EOF
chr1${TAB}5${TAB}10${TAB}myseq${TAB}AAACC
EOF

"$meta_executable" \
  --input_bed "$TMPDIR/test.bed" \
  --input_fasta "$TMPDIR/test.fa" \
  --bed_out \
  --output "$TMPDIR/output.bed"


cmp --silent "$TMPDIR/expected.bed" "$TMPDIR/output.bed" || { echo "Files when using --bed_out are different."; exit 1; }

# Create dummy bed file for strandedness
cat > "$TMPDIR/test_strandedness.bed" <<EOF
chr1${TAB}20${TAB}25${TAB}forward${TAB}1${TAB}+
chr1${TAB}20${TAB}25${TAB}reverse${TAB}1${TAB}-
EOF

# Create expected tab-delimited file for --bed_out
cat > "$TMPDIR/expected_strandedness.fasta" <<EOF
>forward(+)
CGCTA
>reverse(-)
TAGCG
EOF

"$meta_executable" \
  --input_bed "$TMPDIR/test_strandedness.bed" \
  --input_fasta "$TMPDIR/test.fa" \
  -s \
  --name_only \
  --output "$TMPDIR/output_strandedness.fasta"


cmp --silent "$TMPDIR/expected_strandedness.fasta" "$TMPDIR/output_strandedness.fasta" || { echo "Files when using -s are different."; exit 1; }

