#!/bin/bash

# GATK4-specific test helper functions for biobox components
#
# Source this file (in addition to the global /src/_utils/test_helpers.sh)
# via:
#   source "$meta_resources_dir/gatk4/test_helpers.sh"

# Create a samtools-compatible .fai index for a FASTA file without samtools.
#
# Usage: create_test_fasta_fai "/path/to/input.fasta" "/path/to/output.fai"
create_test_fasta_fai() {
  local fasta_path="$1"
  local fai_path="$2"

  log "Creating FASTA index (.fai) for: $fasta_path"

  # Re-implements the samtools .fai algorithm: for each sequence, walk the
  # file byte-by-byte (via line lengths) tracking where the sequence's bases
  # start (offset) and how the bases are wrapped (linebases/linewidth), then
  # emit "name\tseqlen\toffset\tlinebases\tlinewidth" once the next header
  # (or EOF) closes out the record.
  awk '
    BEGIN { name = ""; seqlen = 0; offset = 0; linebases = 0; linewidth = 0; bytepos = 0 }
    {
      # +1 accounts for the newline awk strips from $0
      linelen = length($0) + 1
      if (substr($0, 1, 1) == ">") {
        # Header line: flush the previous record (if any), then start a new
        # one. The sequence name is everything up to the first space/tab.
        if (name != "") {
          print name"\t"seqlen"\t"offset"\t"linebases"\t"linewidth
        }
        name = substr($0, 2)
        sub(/[ \t].*/, "", name)
        seqlen = 0
        linebases = 0
        linewidth = 0
        bytepos += linelen
        offset = bytepos
      } else {
        # Sequence line: linebases/linewidth are fixed from the first
        # sequence line of the record (samtools assumes uniform wrapping),
        # seqlen accumulates the total base count across all its lines.
        if (linebases == 0) {
          linebases = length($0)
          linewidth = linelen
        }
        seqlen += length($0)
        bytepos += linelen
      }
    }
    END {
      # Flush the final record, since there is no following header line to
      # trigger it.
      if (name != "") {
        print name"\t"seqlen"\t"offset"\t"linebases"\t"linewidth
      }
    }
  ' "$fasta_path" > "$fai_path"

  log "✓ Created FASTA index: $fai_path"
}

# Sort a SAM/BAM file and build its .bai index using GATK's bundled Picard
# tools (SortSam + BuildBamIndex)
#
# Usage: sort_and_index_bam "/path/to/input.sam_or_bam" "/path/to/output.bam" [sort_order]
sort_and_index_bam() {
  local input_path="$1"
  local output_bam="$2"
  local sort_order="${3:-coordinate}"

  log "Sorting $input_path -> $output_bam (order: $sort_order)"
  gatk SortSam \
    --INPUT "$input_path" \
    --OUTPUT "$output_bam" \
    --SORT_ORDER "$sort_order" \
    --VERBOSITY ERROR

  log "Building BAM index for: $output_bam"
  gatk BuildBamIndex \
    --INPUT "$output_bam" \
    --OUTPUT "${output_bam%.bam}.bai" \
    --VERBOSITY ERROR

  log "✓ Sorted and indexed BAM: $output_bam (index: ${output_bam%.bam}.bai)"
}

# Create a synthetic reference FASTA together with its .fai and .dict
# companion files (the "reference trio" most GATK4 tools require to find
# indices/dictionaries by basename next to the FASTA).
#
# Usage: create_test_reference "/path/to/reference.fasta" [num_seqs=1] [seq_length=2000]
# Produces "$fasta_path", "${fasta_path}.fai", and a sibling ".dict" file
create_test_reference() {
  local fasta_path="$1"
  local num_seqs="${2:-1}"
  local seq_length="${3:-2000}"
  local dict_path="${fasta_path%.fasta}.dict"

  create_test_fasta "$fasta_path" "$num_seqs" "$seq_length"

  log "Creating sequence dictionary for: $fasta_path"
  gatk CreateSequenceDictionary -R "$fasta_path" -O "$dict_path" --VERBOSITY ERROR

  create_test_fasta_fai "$fasta_path" "${fasta_path}.fai"
}

# Write a SAM file with one read-group-tagged 100bp read (ATCG-repeat
# content, matching the global create_test_fasta helper's output) at each
# given 1-based position. Positions may repeat, e.g. to create duplicate
# reads for a MarkDuplicates test.
#
# Usage: create_test_sam "/path/to/reads.sam" contig_name contig_length \
#          read_group_id sample_name pos1 [pos2 ...]
create_test_sam() {
  local sam_path="$1"
  local contig_name="$2"
  local contig_length="$3"
  local read_group_id="$4"
  local sample_name="$5"
  shift 5
  local positions=("$@")

  local read_length=100
  local read_seq
  read_seq=$(printf 'ATCG%.0s' $(seq 1 $((read_length / 4))))
  local qual
  qual=$(printf 'I%.0s' $(seq 1 "$read_length"))

  log "Writing SAM reads for $sample_name (${#positions[@]} reads) to: $sam_path"

  {
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:%s\tLN:%d\n' "$contig_name" "$contig_length"
    printf '@RG\tID:%s\tSM:%s\tLB:lib1\tPL:ILLUMINA\n' "$read_group_id" "$sample_name"
    local i=0
    local pos
    for pos in "${positions[@]}"; do
      i=$((i + 1))
      printf 'read%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:%s\n' \
        "$i" "$contig_name" "$pos" "$read_length" "$read_seq" "$qual" "$read_group_id"
    done
  } > "$sam_path"

  log "✓ Wrote SAM reads: $sam_path"
}

# Convenience wrapper around create_test_sam: N 100bp reads tiled across the
# contig at a fixed step (default step 50, i.e. consecutive reads overlap by
# half their length).
#
# Usage: create_test_sam_reads "/path/to/reads.sam" contig_name contig_length \
#          read_group_id sample_name [num_reads=6] [step=50]
create_test_sam_reads() {
  local sam_path="$1"
  local contig_name="$2"
  local contig_length="$3"
  local read_group_id="$4"
  local sample_name="$5"
  local num_reads="${6:-6}"
  local step="${7:-50}"

  local positions=()
  local i
  for ((i = 0; i < num_reads; i++)); do
    positions+=("$((i * step + 1))")
  done

  create_test_sam "$sam_path" "$contig_name" "$contig_length" "$read_group_id" "$sample_name" "${positions[@]}"
}

# Create a minimal single-record known-sites VCF at a 1-based position within
# a contig whose sequence is the "ATCG" repeat produced by create_test_fasta,
# and index it with `gatk IndexFeatureFile`. REF/ALT alleles are derived
# automatically from the position using the ATCG cycle, unless explicitly
# overridden.
#
# Usage: create_test_known_sites_vcf "/path/to/known_sites.vcf" contig_name \
#          contig_length position [ref_base] [alt_base]
create_test_known_sites_vcf() {
  local vcf_path="$1"
  local contig_name="$2"
  local contig_length="$3"
  local position="$4"
  local cycle="ATCG"
  local idx=$(( (position - 1) % 4 ))
  local ref_base="${5:-${cycle:$idx:1}}"
  local alt_base="${6:-${cycle:$(( (idx + 1) % 4 )):1}}"

  log "Writing known-sites VCF at ${contig_name}:${position} (${ref_base}>${alt_base}) to: $vcf_path"

  {
    echo "##fileformat=VCFv4.2"
    echo "##contig=<ID=${contig_name},length=${contig_length}>"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    printf '%s\t%d\t.\t%s\t%s\t50\tPASS\t.\n' "$contig_name" "$position" "$ref_base" "$alt_base"
  } > "$vcf_path"

  gatk IndexFeatureFile --input "$vcf_path" --verbosity ERROR

  log "✓ Wrote and indexed known-sites VCF: $vcf_path"
}
