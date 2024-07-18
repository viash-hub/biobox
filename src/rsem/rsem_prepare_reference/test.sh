
#!/bin/bash

set -e pipefail

echo ">>> Testing $meta_functionality_name"

cat > genome.fasta <<'EOF'
>Sheila
GCTAGCTCAGAAAAAAAAAA
EOF

echo ">>> Prepare RSEM reference without gene annotations"
"$meta_executable" \
  --reference_fasta_files genome.fasta \
  --reference_name test \
  --output RSEM_index \

echo ">>> Checking whether output files exist"
[ ! -d "RSEM_index" ] && echo "RSEM index does not exist!" && exit 1
[ ! -f "RSEM_index/test.grp" ] && echo "test.grp does not exist!" && exit 1
[ ! -f "RSEM_index/test.n2g.idx.fa" ] && echo "test.n2g.idx.fa does not exist!" && exit 1
[ ! -f "RSEM_index/test.ti" ] && echo "test.ti does not exist!" && exit 1
[ ! -f "RSEM_index/test.idx.fa" ] && echo "test.idx.fa does not exist!" && exit 1
[ ! -f "RSEM_index/test.seq" ] && echo "test.seq does not exist!" && exit 1
[ ! -f "RSEM_index/test.transcripts.fa" ] && echo "test.transcripts.fa does not exist!" && exit 1

# echo ">>> Checking whether output is correct"
# [ ! -f "RSEM_index/test.grp" ] && echo "test.grp does not exist!" && exit 1
# [ ! -f "RSEM_index/test.n2g.idx.fa" ] && echo "test.n2g.idx.fa does not exist!" && exit 1
# [ ! -f "RSEM_index/test.ti" ] && echo "test.ti does not exist!" && exit 1
# [ ! -f "RSEM_index/test.idx.fa" ] && echo "test.idx.fa does not exist!" && exit 1
# [ ! -f "RSEM_index/test.seq" ] && echo "test.seq does not exist!" && exit 1
# [ ! -f "RSEM_index/test.transcripts.fa" ] && echo "test.transcripts.fa does not exist!" && exit 1