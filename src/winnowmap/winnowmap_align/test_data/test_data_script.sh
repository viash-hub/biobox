#!/bin/bash
# Generate synthetic test data for winnowmap_align.
#
# Run this script once to produce:
#   reference.fasta - 2 contigs: one unique, one carrying tandem-repeat arrays
#   reads.fastq     - 16 simulated long reads (2000 bp) sampled from the reference
#
# The reference deliberately contains two tandem arrays with DIFFERENT copy
# numbers. Winnowmap's repetitive k-mer step is
#   meryl print greater-than distinct=0.9998
# which keeps only the top 0.02% of distinct k-mers by count. A uniformly random
# reference yields an empty set, and so does a single tandem array (all of its
# k-mers share one count, so the quantile threshold lands exactly on them and
# "greater-than" excludes every one). Two arrays of differing copy number put
# the threshold at the lower count and let the higher-count k-mers through, so
# the weighted-minimizer path is actually exercised. This mirrors test.sh.
#
# Usage (from repo root):
#   bash src/winnowmap/winnowmap_align/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ref_file="${SCRIPT_DIR}/reference.fasta"
fastq_file="${SCRIPT_DIR}/reads.fastq"

# ---------------------------------------------------------------------------
# Generate reference.fasta and reads.fastq with Python for reliability
# (bash loops over /dev/urandom are slow and SIGPIPE-prone)
# ---------------------------------------------------------------------------
python3 - "$ref_file" "$fastq_file" <<'PYEOF'
import random
import sys

ref_path, fastq_path = sys.argv[1], sys.argv[2]

random.seed(42)
BASES = "ACGT"

def rand_seq(n):
    return "".join(random.choices(BASES, k=n))

def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

contigs = {}

# Contig 1: unique background sequence
contigs["chr_unique"] = rand_seq(120000)

# Contig 2: unique flanks around two tandem arrays of differing copy number.
# motif_a is high-copy and survives the 0.9998 cutoff; motif_b sets the cutoff.
motif_a = rand_seq(20)
motif_b = rand_seq(31)
contigs["chr_repeat"] = (
    rand_seq(20000)
    + motif_a * 600
    + rand_seq(10000)
    + motif_b * 250
    + rand_seq(20000)
)

with open(ref_path, "w") as f:
    for name, seq in contigs.items():
        f.write(f">{name}\n{wrap(seq)}\n")
        print(f"  {name}: {len(seq)} bp")

# Reads are sampled per contig - never across a contig boundary, which would
# produce chimeric reads that cannot align.
READ_LEN = 2000
READS_PER_CONTIG = 8
SUB_RATE = 0.01

def mutate(seq, rate):
    out = []
    for base in seq:
        out.append(random.choice(BASES) if random.random() < rate else base)
    return "".join(out)

n = 0
with open(fastq_path, "w") as f:
    for name, seq in contigs.items():
        for _ in range(READS_PER_CONTIG):
            start = random.randint(0, len(seq) - READ_LEN)
            read = mutate(seq[start:start + READ_LEN], SUB_RATE)
            n += 1
            f.write(f"@read{n}_{name}_{start}\n{read}\n+\n{'I' * READ_LEN}\n")
print(f"  {n} reads x {READ_LEN} bp written")
PYEOF

echo ""
echo "Done. Run the component with:"
echo ""
echo "  viash run src/winnowmap/winnowmap_align/config.vsh.yaml -- \\"
echo "    --reference ${ref_file} \\"
echo "    --query     ${fastq_file} \\"
echo "    --preset    map-ont \\"
echo "    --output    ${SCRIPT_DIR}/alignment.sam"
echo ""
echo "Or as a coordinate-sorted, indexed BAM:"
echo ""
echo "  viash run src/winnowmap/winnowmap_align/config.vsh.yaml -- \\"
echo "    --reference    ${ref_file} \\"
echo "    --query        ${fastq_file} \\"
echo "    --preset       map-ont \\"
echo "    --bam \\"
echo "    --output       ${SCRIPT_DIR}/alignment.bam \\"
echo "    --output_index ${SCRIPT_DIR}/alignment.bam.bai"
