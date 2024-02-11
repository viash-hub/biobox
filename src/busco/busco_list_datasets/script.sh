#!/bin/bash

## VIASH START
## VIASH END

busco --list-datasets | awk '/^#{40}/{flag=1; next} flag{print}' > $par_output