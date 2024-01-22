#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run arriba with blacklist"
"$meta_executable" \
  ...

echo ">> Checking output"
...

echo "> Test successful"