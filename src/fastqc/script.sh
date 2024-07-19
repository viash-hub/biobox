#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_casava" == "false" ]] && unset par_

# run fastqc
fastqc \
  ${par_casava:+--casava} \
  $par_input
  