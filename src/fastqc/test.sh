#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/fastqc"
meta_resources_dir="src/fastqc"
## VIASH END

# Create directories for tests
echo "Creating Test Data..."
mkdir test_data

# Create and populate input.fasta
cat > "test_data/input.fq" <<EOL
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
EOL

# Test 1: Run fastqc with default parameters


