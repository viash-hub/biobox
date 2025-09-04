#!/bin/bash

TOOL=bcftools
DOCKER_IMAGE="quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

SUBCOMMANDS=(
  # indexing
  index

  # vcf/bcf manipulation
  annotate
  concat
  convert
  head
  isec
  merge
  norm
  plugin
  query
  reheader
  sort
  view

  # vcf/bcf analysis
  call
  consensus
  cnv
  csq
  filter
  gtcheck
  mpileup
  polysomy
  roh
  stats
)

for SUBCOMMAND in "${SUBCOMMANDS[@]}"; do

  DIR="src/$TOOL/${TOOL}_$SUBCOMMAND"
  DEST="$DIR/help.txt"
  CFG="$DIR/config.vsh.yaml"
  CMD="docker run --rm $DOCKER_IMAGE $TOOL $SUBCOMMAND --help 2>&1 | grep -v unrecognized"

  # if config.vsh.yaml does not exist in dir, skip
  if [ ! -f "$CFG" ]; then
    echo "Config file $CFG does not exist, skipping."
    continue
  fi
  
  echo "Generating help for $TOOL $SUBCOMMAND"

  # # create dir if not exists
  # mkdir -p "$(dirname "$DEST")"

  # add header
  printf '```bash\n%s\n```\n' "$CMD" > "$DEST"

  # add help to file
  eval "$CMD" >> "$DEST" 2>&1
done
