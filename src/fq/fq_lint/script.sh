#!/bin/bash

## VIASH START
# par_input_1="reads_1.fastq.gz"
# par_input_2="reads_2.fastq.gz"
# par_lint_mode="panic"
# par_disable_validator="S001;P002"
## VIASH END

# Exit immediately if a command exits with a non-zero status.
set -eo pipefail

# split the disable_validator string into an array
IFS=';' read -r -a par_disable_validator <<< "$par_disable_validator"

# Construct and execute the fq lint command.
fq lint \
  ${par_lint_mode:+--lint-mode "$par_lint_mode"} \
  ${par_single_read_validation_level:+--single-read-validation-level "$par_single_read_validation_level"} \
  ${par_paired_read_validation_level:+--paired-read-validation-level "$par_paired_read_validation_level"} \
  ${par_record_definition_separator:+--record-definition-separator "$par_record_definition_separator"} \
  ${par_disable_validator[@]/#/--disable-validator } \
  "$par_input_1" \
  ${par_input_2:+"$par_input_2"}
