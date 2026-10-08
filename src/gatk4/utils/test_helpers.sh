#!/bin/bash

# GATK4-specific test helper functions for biobox components
#
# Source this file (in addition to the global /src/_utils/test_helpers.sh)
# via:
#   source "$meta_resources_dir/gatk4/test_helpers.sh"

# Create a .fai index for a FASTA file using samtools.
#
# Usage: create_test_fasta_fai "/path/to/input.fasta" "/path/to/output.fai"
create_test_fasta_fai() {
  local fasta_path="$1"
  local fai_path="$2"

  log "Creating FASTA index (.fai) for: $fasta_path"

  samtools faidx "$fasta_path"

  if [[ "$fai_path" != "${fasta_path}.fai" ]]; then
    mv "${fasta_path}.fai" "$fai_path"
  fi

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

# Count the number of reads in a SAM/BAM/CRAM file using GATK's CountReads
# tool.
#
# Usage: count_bam_reads "/path/to/file.bam"
count_bam_reads() {
  local bam_path="$1"
  gatk CountReads --input "$bam_path" --verbosity ERROR 2>/dev/null | tail -n 1
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

# Somatic variant calling fixtures
#
# A 500bp single-contig ("chr1") reference with a known SNV at 1-based
# position 250 (REF T > ALT A). Read pairs are two 100bp mates at 1-based
# positions 200 and 350. The pos-200 mate spans the SNV (it is the 51st base
# of that mate), the pos-350 mate never carries the variant.
TEST_SOMATIC_SEQ="AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATAAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAAACAGCTCAGGAGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATACGTAGGAGGTGTGCCTTCGTAGGCTGTTTCTCAGGACGCCCAACTATTCTTTCCAATCCTACATCTGTTTCTTGCGTCGTAGCGGGACCCTCCATTGTTACTTATTAGGTTCTCGTTATGTCTCATAATCTC"
TEST_SOMATIC_SEQ_REF="AAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTA"
TEST_SOMATIC_SEQ_ALT="AAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGACGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTA"
TEST_SOMATIC_SEQ_MATE="TACGACGCAAGAAACAGATGTAGGATTGGAAAGAATAGTTGGGCGTCCTGAGAAACAGCCTACGAAGGCACACCTCCTACGTATCGATACCATTCGCTTC"

# Create the 500bp somatic test reference FASTA together with its .fai and
# .dict companion files.
#
# Usage: create_test_somatic_reference "/path/to/reference.fasta"
# Produces "$fasta_path", "${fasta_path}.fai", and a sibling ".dict" file
create_test_somatic_reference() {
  local fasta_path="$1"
  local dict_path="${fasta_path%.fasta}.dict"

  log "Writing somatic test reference to: $fasta_path"
  printf '>chr1\n%s\n' "$TEST_SOMATIC_SEQ" > "$fasta_path"

  log "Creating sequence dictionary for: $fasta_path"
  gatk CreateSequenceDictionary -R "$fasta_path" -O "$dict_path" --VERBOSITY ERROR

  create_test_fasta_fai "$fasta_path" "${fasta_path}.fai"
}

# Create a sorted and indexed BAM of 8 read pairs over the somatic test SNV.
# With mode "ref" all reads carry the reference allele (a matched normal
# sample). With mode "mixed" the reads alternate between the reference and
# alternate alleles (a tumor sample with a ~50% allele fraction).
#
# Usage: create_test_somatic_bam "/path/to/output.bam" sample_name read_group_id ref|mixed
# Produces "$bam_path" and "${bam_path%.bam}.bai"
create_test_somatic_bam() {
  local bam_path="$1"
  local sample_name="$2"
  local read_group_id="$3"
  local mode="$4"
  local sam_path="${bam_path%.bam}.sam"
  local qual
  qual=$(printf 'I%.0s' $(seq 1 100))

  log "Writing somatic test reads for $sample_name (mode: $mode) to: $sam_path"
  {
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:chr1\tLN:500\n'
    printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tLB:lib_%s\n' "$read_group_id" "$sample_name" "$read_group_id"
    local i variant_seq
    for i in $(seq 0 7); do
      if [[ "$mode" == "mixed" ]] && ((i % 2 == 0)); then
        variant_seq="$TEST_SOMATIC_SEQ_ALT"
      else
        variant_seq="$TEST_SOMATIC_SEQ_REF"
      fi
      # Alternate which mate is the forward-strand read
      if ((i % 4 < 2)); then
        printf 'pair%d\t99\tchr1\t200\t60\t100M\t=\t350\t250\t%s\t%s\tRG:Z:%s\n' "$i" "$variant_seq" "$qual" "$read_group_id"
        printf 'pair%d\t147\tchr1\t350\t60\t100M\t=\t200\t-250\t%s\t%s\tRG:Z:%s\n' "$i" "$TEST_SOMATIC_SEQ_MATE" "$qual" "$read_group_id"
      else
        printf 'pair%d\t83\tchr1\t350\t60\t100M\t=\t200\t-250\t%s\t%s\tRG:Z:%s\n' "$i" "$TEST_SOMATIC_SEQ_MATE" "$qual" "$read_group_id"
        printf 'pair%d\t163\tchr1\t200\t60\t100M\t=\t350\t250\t%s\t%s\tRG:Z:%s\n' "$i" "$variant_seq" "$qual" "$read_group_id"
      fi
    done
  } > "$sam_path"

  sort_and_index_bam "$sam_path" "$bam_path"
}

# Create the full set of somatic fixtures in a directory:
#
# - the reference (reference.fasta/.fai/.dict)
# - a tumor BAM (tumor.bam/.bai, sample "tumor")
# - a matched normal BAM (normal.bam/.bai, sample "normal").
#
# Then run Mutect2 directly on the tumor/normal pair to produce an unfiltered
# VCF (unfiltered.vcf), its stats file (unfiltered.vcf.stats) and an F1R2
# tarball (f1r2.tar.gz) for testing the downstream somatic tools.
#
# Usage: create_test_mutect2_outputs "/path/to/dir"
create_test_mutect2_outputs() {
  local out_dir="$1"

  create_test_somatic_reference "$out_dir/reference.fasta"
  create_test_somatic_bam "$out_dir/tumor.bam" tumor rg_tumor mixed
  create_test_somatic_bam "$out_dir/normal.bam" normal rg_normal ref

  log "Running Mutect2 on the somatic test tumor/normal pair..."
  gatk Mutect2 \
    --reference "$out_dir/reference.fasta" \
    --input "$out_dir/tumor.bam" \
    --input "$out_dir/normal.bam" \
    --normal-sample normal \
    --output "$out_dir/unfiltered.vcf" \
    --f1r2-tar-gz "$out_dir/f1r2.tar.gz" \
    --verbosity ERROR

  log "✓ Created Mutect2 outputs in: $out_dir"
}
