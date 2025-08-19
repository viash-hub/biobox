#!/bin/bash

TOOL=bedtools
DOCKER_IMAGE="quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"

SUBCOMMANDS=(
  # genome arithmetic
  intersect
  window
  closest
  coverage
  map
  genomecov
  merge
  cluster
  complement
  shift
  subtract
  slop
  flank
  sort
  random
  shuffle
  sample
  spacing
  annotate
  # multi-way file comparisons
  multiinter
  unionbedg
  # paired-end manipulation
  pairtobed
  pairtopair
  # format conversion
  bamtobed
  bedtobam
  bamtofastq
  bedpetobam
  # fasta manipulation tools
  bed12tobed6
  getfasta
  maskfasta
  # bam focused tools
  multicov
  tag
  # statistical relationships
  jaccard
  reldist
  fisher
  # miscellaneous tools
  overlap
  igv
  links
  makewindows
  groupby
  expand
  split
  summary
)

for SUBCOMMAND in "${SUBCOMMANDS[@]}"; do
  echo "Generating help for $TOOL $SUBCOMMAND"

  DEST="src/$TOOL/${TOOL}_$SUBCOMMAND/help.txt"
  CMD="docker run --rm \"$DOCKER_IMAGE\" \"$TOOL\" \"$SUBCOMMAND\" -h"

  # skip if help doesn't already exist
  if [ ! -f "$DEST" ]; then
    continue
  fi

  # create dir if not exists
  mkdir -p "$(dirname "$DEST")"

  # add header
  printf '```bash\n%s\n```\n' "$CMD" > "$DEST"

  # add help to file
  eval "$CMD" >> "$DEST" 2>&1
done
