#!/bin/bash

## VIASH START
## VIASH END


# unset flags
[[ "$par_emblmygff3" == "false" ]] && unset par_emblmygff3
[[ "$par_discard" == "false" ]] && unset par_discard
[[ "$par_keep" == "false" ]] && unset par_keep

# replace ';' with ','
par_primary_tag=$(echo $par_primary_tag | tr ';' ',')

# run agat_convert_embl2gff
agat_convert_embl2gff.pl \
  --embl "$par_embl" \
  -o "$par_output" \
  ${par_emblmygff3:+--emblmygff3} \
  ${par_primary_tag:+--primary_tag "${par_primary_tag}"} \
  ${par_discard:+-d} \
  ${par_keep:+-k} \
  ${par_config:+--config "${par_config}"}
