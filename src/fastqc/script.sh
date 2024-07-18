#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_" == "false" ]] && unset par_

# run fastqc
fastqc \
  ${par_:+ "${par_}"} \
  $par_input
  