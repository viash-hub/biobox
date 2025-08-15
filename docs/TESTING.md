# Testing Guide

This guide covers best practices for writing comprehensive test scripts for biobox components.

> **📌 Important:** All new test scripts should use the **centralized test helpers** located at `src/_utils/test_helpers.sh`. This eliminates code duplication and ensures consistency across all components.

## Table of Contents

- [Core Principles](#core-principles)
- [Test Script Structure](#test-script-structure)
- [Centralized Test Helpers](#centralized-test-helpers)
- [Test Scenarios](#test-scenarios)
- [Best Practices](#best-practices)
- [Viash Testing Features](#viash-testing-features)
- [Static Test Data](#static-test-data)

## Core Principles

### 1. Generate Test Data in Scripts

**Preferred approach:** Generate test data within the test script using the centralized helper functions.

```bash
# Generate test data using centralized helpers
create_test_fasta "$meta_temp_dir/input.fasta" 3 50
create_test_fastq "$meta_temp_dir/reads.fastq" 10 35
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
  - type: file
    path: /src/_utils/test_helpers.sh
```

Only add static test files if absolutely necessary:

```yaml
test_resources:
  - type: bash_script
    path: test.sh
  - type: file
    path: /src/_utils/test_helpers.sh
  - type: file
    path: test_data  # Only if data generation is impractical
```

## Test Script Structure

### Configuration Setup

Add the test helpers as a resource in your component configuration:

```yaml
test_resources:
  - type: bash_script
    path: test.sh
  - type: file
    path: /src/_utils/test_helpers.sh
```

### Basic Test Template

```bash
#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# --- Test Case 1: Basic functionality ---
log "Starting TEST 1: Basic functionality"

# Create and validate test data
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"
create_test_fasta "$test_data_dir/input.fasta" 3 50
check_file_exists "$test_data_dir/input.fasta" "input FASTA file"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/input.fasta" \
  --output "$meta_temp_dir/test1"

log "Validating TEST 1 outputs..."
check_dir_exists "$meta_temp_dir/test1" "output directory"
check_file_exists "$meta_temp_dir/test1/result.txt" "result file"
check_file_not_empty "$meta_temp_dir/test1/result.txt" "result file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Advanced parameters ---
log "Starting TEST 2: Advanced parameters"

# Create different test data
create_test_fastq "$test_data_dir/input.fastq" 10 35
check_file_exists "$test_data_dir/input.fastq" "input FASTQ file"

log "Executing $meta_name with advanced parameters..."
"$meta_executable" \
  --input "$test_data_dir/input.fastq" \
  --output "$meta_temp_dir/test2" \
  --threads 2 \
  --verbose

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/test2/advanced_result.txt" "advanced result file"
check_file_contains "$meta_temp_dir/test2/advanced_result.txt" "expected_pattern" "advanced result file"

log "✅ TEST 2 completed successfully"

print_test_summary "All tests completed successfully"
```

## Centralized Test Helpers

The centralized test helpers located at `src/_utils/test_helpers.sh` provide comprehensive testing functionality to ensure consistency across all biobox components.

### Available Functions

#### Logging Functions
- `log "message"` - Log with timestamp
- `log_warn "message"` - Warning message  
- `log_error "message"` - Error message

#### File/Directory Validation
- `check_file_exists path "description"` - Verify file exists
- `check_dir_exists path "description"` - Verify directory exists
- `check_file_not_exists path "description"` - Verify file doesn't exist
- `check_dir_not_exists path "description"` - Verify directory doesn't exist
- `check_file_empty path "description"` - Verify file is empty
- `check_file_not_empty path "description"` - Verify file is not empty

#### Content Validation
- `check_file_contains path "text" "description"` - Verify file contains text
- `check_file_not_contains path "text" "description"` - Verify file doesn't contain text
- `check_file_matches_regex path "pattern" "description"` - Verify file matches regex
- `check_file_line_count path count "description"` - Verify line count

#### Test Data Generation
- `create_test_fasta path [num_seqs] [seq_length]` - Generate FASTA file
- `create_test_fastq path [num_reads] [read_length]` - Generate FASTQ file
- `create_test_gtf path [num_genes]` - Generate GTF file
- `create_test_gff path [num_features]` - Generate GFF file
- `create_test_bed path [num_intervals]` - Generate BED file
- `create_test_csv path [num_rows]` - Generate CSV file
- `create_test_tsv path [num_rows]` - Generate TSV file

#### Utility Functions
- `setup_test_env` - Initialize test environment with strict error handling
- `print_test_summary "test_name"` - Print completion message

### Usage Example

```bash
#!/bin/bash

## VIASH START
## VIASH END

# Source centralized helpers
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

# Generate test data
create_test_fasta "$meta_temp_dir/input.fasta" 3 50
check_file_exists "$meta_temp_dir/input.fasta" "input FASTA file"

# Run component
"$meta_executable" \
  --input "$meta_temp_dir/input.fasta" \
  --output "$meta_temp_dir/output.txt"

# Validate output
check_file_exists "$meta_temp_dir/output.txt" "result file"
check_file_contains "$meta_temp_dir/output.txt" "expected_pattern" "result file"

print_test_summary "Basic functionality test"
```

## Test Scenarios

### 1. Basic Functionality

Test the component with minimal, essential parameters:

```bash
log "Starting TEST 1: Basic functionality"

create_test_fasta "$meta_temp_dir/input.fasta" 3 50

"$meta_executable" \
  --input "$meta_temp_dir/input.fasta" \
  --output "$meta_temp_dir/output.txt"

check_file_exists "$meta_temp_dir/output.txt" "output file"
check_file_not_empty "$meta_temp_dir/output.txt" "output file"

log "✅ TEST 1 completed successfully"
```

### 2. Multiple Input Files

Test with multiple input files or complex input scenarios:

```bash
log "Starting TEST 2: Multiple input files"

create_test_fasta "$meta_temp_dir/input1.fasta" 2 30
create_test_fasta "$meta_temp_dir/input2.fasta" 2 30

"$meta_executable" \
  --input "$meta_temp_dir/input1.fasta;$meta_temp_dir/input2.fasta" \
  --output "$meta_temp_dir/output.txt"

check_file_exists "$meta_temp_dir/output.txt" "merged output file"

log "✅ TEST 2 completed successfully"
```

### 3. Optional Parameters

Test with optional parameters and advanced features:

```bash
log "Starting TEST 3: Optional parameters"

create_test_fastq "$meta_temp_dir/input.fastq" 10 35

"$meta_executable" \
  --input "$meta_temp_dir/input.fastq" \
  --output "$meta_temp_dir/output.txt" \
  --threads 2 \
  --verbose

check_file_exists "$meta_temp_dir/output.txt" "output file with options"
check_file_contains "$meta_temp_dir/output.txt" "verbose" "verbose output"

log "✅ TEST 3 completed successfully"
```

### 4. Edge Cases

Test with edge cases like empty files or unusual inputs:

```bash
log "Starting TEST 4: Edge case - empty input"

# Create empty input file
touch "$meta_temp_dir/empty.fasta"

# Test should handle empty input gracefully
if "$meta_executable" \
  --input "$meta_temp_dir/empty.fasta" \
  --output "$meta_temp_dir/output.txt" 2>/dev/null; then
  log_warn "Component succeeded with empty input - checking output"
  check_file_exists "$meta_temp_dir/output.txt" "output file for empty input"
else
  log "Expected behavior: Component properly rejected empty input"
fi

log "✅ TEST 4 completed successfully"
```

### 5. Error Handling

Test proper error handling for invalid inputs:

```bash
log "Starting TEST 5: Error handling"

# Test with non-existent input file
if "$meta_executable" \
  --input "/non/existent/file.txt" \
  --output "$meta_temp_dir/output.txt" 2>/dev/null; then
  log_error "Component should have failed with non-existent input"
  exit 1
else
  log "✅ Component properly handled non-existent input file"
fi

log "✅ TEST 5 completed successfully"
```

## Best Practices

### 1. Use Centralized Test Helpers

Always use the centralized test helpers instead of defining functions individually:

```bash
# ✅ Recommended: Use centralized helpers
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

# ❌ NOT recommended: Defining functions individually
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') [TEST] $*"; }
```

### 2. Strict Error Handling

The centralized helpers automatically provide strict error handling via `setup_test_env`:

```bash
# Automatically enabled by setup_test_env:
set -euo pipefail  # Exit on errors, undefined variables, pipe failures
export LC_ALL=C    # Consistent locale for reproducible results
```

### 3. Descriptive Validation

Use descriptive validation functions with meaningful descriptions:

```bash
# ✅ Good: Descriptive validation
check_file_exists "$output_file" "filtered feature matrix"
check_file_not_exists "$bam_file" "BAM file (should be disabled by default)"
check_file_contains "$result_file" "expected_pattern" "analysis results"

# ❌ Less helpful: Basic validation without context
check_file_exists "$output_file"
```

### 4. Organized Structure

Use `$meta_temp_dir` and create organized test structure:

```bash
# Create organized test structure
test_data_dir="$meta_temp_dir/test_data"
test_output_dir="$meta_temp_dir/test_output"
mkdir -p "$test_data_dir" "$test_output_dir"

create_test_fasta "$test_data_dir/input.fasta" 3 50
```

### 5. Clear Test Output

Use consistent logging with clear test boundaries:

```bash
log "Starting TEST 1: Basic functionality"
log "Executing $meta_name with basic parameters..."
log "Validating TEST 1 outputs..."
log "✅ TEST 1 completed successfully"

# Final summary
print_test_summary "All tests completed successfully"
```

### 6. Comprehensive Content Validation

Don't just check that files exist - validate their content:

```bash
# Check existence and content
check_file_exists "$meta_temp_dir/output.txt" "analysis results"
check_file_not_empty "$meta_temp_dir/output.txt" "analysis results"
check_file_contains "$meta_temp_dir/output.txt" "Number of sequences" "result summary"
check_file_line_count "$meta_temp_dir/output.txt" 10 "expected number of results"
```

### 7. Multiple Test Scenarios

Include comprehensive test coverage:

```bash
# Test 1: Basic functionality
log "Starting TEST 1: Basic functionality"
# ... test implementation ...
log "✅ TEST 1 completed successfully"

# Test 2: Advanced options
log "Starting TEST 2: Advanced options"
# ... test implementation ...
log "✅ TEST 2 completed successfully"

# Test 3: Edge cases
log "Starting TEST 3: Edge case handling"
# ... test implementation ...
log "✅ TEST 3 completed successfully"

print_test_summary "All tests completed successfully"
```

## Viash Testing Features

### Running Tests

```bash
# Test a single component
viash test config.vsh.yaml

# Test with specific resources
viash test config.vsh.yaml --cpus 4 --memory 8GB

# Test with specific setup strategy
viash test config.vsh.yaml --setup build --verbose

# Keep temporary files for debugging
viash test config.vsh.yaml --keep true

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

### Meta Variables in Tests

Your test scripts automatically have access to important meta variables:

- `$meta_executable` - Path to the built component executable
- `$meta_temp_dir` - Temporary directory for test files (automatically cleaned up)
- `$meta_name` - Component name for logging
- `$meta_resources_dir` - Path to test resources

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
  - type: file
    path: /src/_utils/test_helpers.sh
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

## Static Test Data

### When to Use Static Test Data

Only use static test files when:

- The tool requires very specific, complex file formats that are difficult to generate
- Generating equivalent test data is impractical or overly complex
- You need real-world data to validate complex algorithms
- Test data is very small (<1KB preferred, <10KB maximum)

### Guidelines for Static Test Data

If you must use static test data:

1. **Keep files small** - Prefer <1KB, maximum <10KB
2. **Document the source** - How was it created?
3. **Use minimal examples** - Strip down to essential features
4. **Consider alternatives** - Can you generate equivalent data?

```bash
# test_data/README.md
# Test data for complex_tool component
# Source: https://github.com/example/dataset
# Generated with: tool --export-sample --format minimal
# Date: 2025-01-01
# Size: 847 bytes
# Purpose: Tests complex file format parsing
```

### Referencing Static Test Data

```yaml
test_resources:
  - type: bash_script
    path: test.sh
  - type: file
    path: /src/_utils/test_helpers.sh
  - type: file
    path: test_data
```

```bash
# In your test script
static_data="$meta_resources_dir/test_data/sample.complex"
check_file_exists "$static_data" "static test data"

"$meta_executable" --input "$static_data" --output "$meta_temp_dir/output.txt"
```
