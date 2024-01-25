#!/bin/bash

cat <<-END_VERSIONS
busco: "$(busco --version 2>&1 | sed -n 's/BUSCO \([0-9.]*\)/\1/p')"
END_VERSIONS