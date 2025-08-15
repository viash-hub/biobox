# Testing Guide

This guide covers best practices for writing comprehensive test scripts for biobox components.

## Table of Contents
- [Core Principles](#core-principles)
- [Test Script Structure](#test-script-structure)
- [Helper Functions](#helper-functions)
- [Test Scenarios](#test-scenarios)
- [Best Practices](#best-practices)
- [Viash Testing Features](#viash-testing-features)
- [When to Use Static Test Data](#when-to-use-static-test-data)

## Core Principles

### 1. Generate Test Data in Scripts

**Preferred approach:**
```bash
# Generate test data within the test script
create_test_fasta() {
  file_path="$1"
  
  cat << 'EOF' > "$file_path"
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
EOF
}
```

**Avoid:**
- Storing static test files in the repository
- Fetching test data from external sources
- Large test datasets

### 2. Self-Contained Tests

Tests should be completely self-contained and not depend on external resources:

```yaml
test_resources:
  - type: bash_script
    path: test.sh
```

Only add static test files if absolutely necessary:

```yaml
test_resources:
  - type: bash_script
    path: test.sh
  - type: file
    path: test_data  # Only if data generation is impractical
```

## Test Script Structure

### Basic Template

```bash
#!/bin/bash

set -e

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_dir_exists() {
  [ -d "$1" ] || { echo "Directory '$1' does not exist" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
#############################################

# --- Helper function to create test data ---
create_test_fasta() {
  file_path="$1"
  
  cat << 'EOF' > "$file_path"
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
EOF
}

# --- Test Case 1: Basic functionality ---
echo ">>> Test 1: Basic functionality"
create_test_fasta "$meta_temp_dir/input.fasta"

echo ">> Running $meta_name..."
"$meta_executable" \
  --input "$meta_temp_dir/input.fasta" \
  --output "$meta_temp_dir/output"

echo ">> Checking output exists..."
assert_dir_exists "$meta_temp_dir/output"

echo ">> Checking output is not empty..."
assert_file_not_empty "$meta_temp_dir/output/result.txt"

echo ">> Checking output content..."
assert_file_contains "$meta_temp_dir/output/result.txt" "expected_pattern"

echo "> All tests succeeded!"
```

## Helper Functions

### Standard Assertions

```bash
# File existence
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}

assert_file_not_exists() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}

# File content
assert_file_empty() {
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}

assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}

# Content matching
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}

assert_file_contains_regex() {
  grep -q -E "$2" "$1" || { echo "File '$1' does not match regex '$2'" && exit 1; }
}

# Directory checks
assert_dir_exists() {
  [ -d "$1" ] || { echo "Directory '$1' does not exist" && exit 1; }
}

# Numeric comparisons
assert_file_line_count() {
  local expected="$2"
  local actual=$(wc -l < "$1")
  [ "$actual" -eq "$expected" ] || { 
    echo "File '$1' has $actual lines, expected $expected" && exit 1; 
  }
}
```

### Data Generation Functions

```bash
# FASTA files
create_test_fasta() {
  local file_path="$1"
  local num_seqs="${2:-2}"
  local seq_length="${3:-64}"
  
  for i in $(seq 1 "$num_seqs"); do
    echo ">seq$i" >> "$file_path"
    # Generate random DNA sequence
    head -c "$seq_length" < /dev/zero | tr '\0' 'A' >> "$file_path"
    echo >> "$file_path"
  done
}

# FASTQ files
create_test_fastq() {
  local file_path="$1"
  local num_reads="${2:-4}"
  
  for i in $(seq 1 "$num_reads"); do
    echo "@read$i" >> "$file_path"
    echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCG" >> "$file_path"
    echo "+" >> "$file_path"
    echo "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> "$file_path"
  done
}

# GTF/GFF files
create_test_gtf() {
  local file_path="$1"
  
  cat << 'EOF' > "$file_path"
chr1	test	gene	1000	2000	.	+	.	gene_id "gene1"
chr1	test	exon	1000	1500	.	+	.	gene_id "gene1"; transcript_id "transcript1"
chr1	test	exon	1600	2000	.	+	.	gene_id "gene1"; transcript_id "transcript1"
EOF
}
```

## Test Scenarios

### 1. Basic Functionality

```bash
echo ">>> Test 1: Basic functionality"
create_test_input "$meta_temp_dir/input.txt"

"$meta_executable" \
  --input "$meta_temp_dir/input.txt" \
  --output "$meta_temp_dir/output.txt"

assert_file_exists "$meta_temp_dir/output.txt"
assert_file_not_empty "$meta_temp_dir/output.txt"
```

### 2. Multiple Input Files

```bash
echo ">>> Test 2: Multiple input files"
create_test_input "$meta_temp_dir/input1.txt"
create_test_input "$meta_temp_dir/input2.txt"

"$meta_executable" \
  --input "$meta_temp_dir/input1.txt;$meta_temp_dir/input2.txt" \
  --output "$meta_temp_dir/output.txt"

assert_file_exists "$meta_temp_dir/output.txt"
```

### 3. Optional Parameters

```bash
echo ">>> Test 3: Optional parameters"
create_test_input "$meta_temp_dir/input.txt"

"$meta_executable" \
  --input "$meta_temp_dir/input.txt" \
  --output "$meta_temp_dir/output.txt" \
  --threads 2 \
  --verbose

assert_file_exists "$meta_temp_dir/output.txt"
```

### 4. Edge Cases

```bash
echo ">>> Test 4: Empty input"
touch "$meta_temp_dir/empty.txt"

"$meta_executable" \
  --input "$meta_temp_dir/empty.txt" \
  --output "$meta_temp_dir/output.txt" || {
  echo "Expected failure with empty input"
}
```

### 5. Error Handling

```bash
echo ">>> Test 5: Error handling - non-existent input"
if "$meta_executable" \
  --input "/non/existent/file.txt" \
  --output "$meta_temp_dir/output.txt" 2>/dev/null; then
  echo "Expected error but command succeeded" && exit 1
fi
echo ">> OK: Properly handled non-existent input file error."
```

## Best Practices

### 1. Use Temporary Directory

Always use `$meta_temp_dir` for temporary files:

```bash
TEMP_DIR="$meta_temp_dir"
create_test_input "$meta_temp_dir/input.txt"
```

### 2. Clear Test Output

Use descriptive echo statements:

```bash
echo ">>> Test 1: Basic functionality"
echo ">> Running $meta_name..."
echo ">> Checking output exists..."
echo ">> OK: Basic functionality test passed."
```

### 3. Test Output Content

Don't just check that files exist - validate their content:

```bash
# Check specific content
assert_file_contains "$meta_temp_dir/output.txt" "expected_result"

# Check file format
assert_file_contains_regex "$meta_temp_dir/output.gtf" "^chr[0-9X-Y]+\t.*\tgene\t"

# Check line count
assert_file_line_count "$meta_temp_dir/output.txt" 10
```

### 4. Test Multiple Scenarios

Include tests for:
- Basic functionality
- Optional parameters
- Multiple inputs
- Edge cases (empty files, etc.)
- Error conditions
- Different output formats

### 5. Fail Fast

Use `set -e` and proper assertions that exit immediately on failure:

```bash
set -e  # Exit on any error

# Each assertion should exit on failure
assert_file_exists "$output_file"
```

## Viash Testing Features

### Running Tests

**Basic test execution:**
```bash
# Test a single component
viash test config.vsh.yaml

# Test with specific resources
viash test config.vsh.yaml --cpus 4 --memory 8GB

# Test with specific setup strategy
viash test config.vsh.yaml --setup build --verbose

# Keep temporary files for debugging
viash test config.vsh.yaml --keep true
```

**Namespace-level testing:**
```bash
# Test all components in parallel
viash ns test --parallel

# Test specific namespace
viash ns test -q alignment --parallel
```

### Test Execution Flow

When running `viash test`, Viash automatically:

1. **Creates temporary directory** (available as `$meta_temp_dir`)
2. **Builds the main executable** 
3. **Builds/pulls Docker image** (if using Docker engine)
4. **Iterates over all test scripts** in `test_resources`
5. **Builds each test into executable** and runs it
6. **Cleans up** temporary files (unless `--keep true`)
7. **Returns exit code 0** if all tests succeed

### Multiple Test Scripts

You can add multiple test scripts to cover different scenarios:

```yaml
test_resources:
  - type: bash_script
    path: test_basic.sh
  - type: bash_script 
    path: test_edge_cases.sh
  - type: bash_script
    path: test_large_data.sh
```

**Important:** All test scripts must pass for the component to be considered tested successfully.

### Testing Different Engines

```bash
# Test with specific engine
viash test config.vsh.yaml --engine docker
viash test config.vsh.yaml --engine native

# Test with specific runner
viash test config.vsh.yaml --runner executable
viash test config.vsh.yaml --runner nextflow
```

### Meta Variables in Tests

Your test scripts automatically have access to important meta variables:

- **`$meta_executable`**: Path to the built component executable
- **`$meta_temp_dir`**: Temporary directory for test files (automatically cleaned up)
- **`$meta_name`**: Component name for logging
- **`$meta_resources_dir`**: Path to test resources

```bash
#!/bin/bash

## VIASH START
## VIASH END

set -e

echo "Testing $meta_name"
echo "Executable: $meta_executable"
echo "Temp dir: $meta_temp_dir"

# Create test data in temp directory
create_test_data "$meta_temp_dir/input.txt"

# Run the component
"$meta_executable" --input "$meta_temp_dir/input.txt" --output "$meta_temp_dir/output.txt"
```

### Advanced Testing Options

```bash
# Test with different container setup strategies
viash test config.vsh.yaml --setup cachedbuild  # Use cached layers (faster)
viash test config.vsh.yaml --setup build        # Clean build from scratch
viash test config.vsh.yaml --setup alwaysbuild  # Always rebuild container

# Test with configuration modifications
viash test config.vsh.yaml -c '.engines[0].image = "ubuntu:22.04"'

# Test with debug mode for troubleshooting
viash test config.vsh.yaml --keep true --verbose
```

For more details, see the [Viash Unit Testing Documentation](https://viash.io/guide/component/unit-testing.html).

## When to Use Static Test Data

Only use static test files when:
- The tool requires very specific, complex file formats
- Generating equivalent test data is impractical
- You need real-world data to validate complex algorithms
- Test data is very small (<1KB preferred, <10KB maximum)

If you must use static test data, document how it was created:

```bash
# test_data/README.md
# Test data source: https://github.com/example/source
# Generated with: command --example --output test.txt
# Date: 2025-01-01
```
