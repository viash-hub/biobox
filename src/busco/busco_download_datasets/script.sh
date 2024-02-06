#!/bin/bash

## VIASH START
## VIASH END


if [ ! -d "$par_download_path" ]; then
    mkdir -p "$par_download_path"
fi

busco \
    --download_path "$par_download_path" \
    --download "$par_download"
    
