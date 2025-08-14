# Script Development Guide

This guide covers best practices for writing runner scripts in biobox components.

## Modern Script Structure

### Basic Template

```bash
#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_option1" == "false" ]] && unset par_option1
[[ "$par_option2" == "false" ]] && unset par_option2

# Build command arguments array
cmd_args=(
    --input "$par_input"
    --output "$par_output"
    ${par_option1:+--option1}
    ${par_option2:+--option2}
    ${par_threads:+--threads "$par_threads"}
    ${par_memory:+--memory "$par_memory"}
)

# Execute command
xxx "${cmd_args[@]}"
```

## Key Principles

### 1. Error Handling

Always use `set -eo pipefail`:
- `set -e`: Exit immediately if a command exits with a non-zero status
- `set -o pipefail`: Exit if any command in a pipeline fails

### 2. Array-Based Arguments

**Preferred approach:**
```bash
cmd_args=(
    --input "$par_input"
    --output "$par_output"
    ${par_option:+--option "$par_option"}
)

xxx "${cmd_args[@]}"
```

**Avoid repetitive appending:**
```bash
# Don't do this
cmd_args+=("--input")
cmd_args+=("$par_input")
cmd_args+=("--output")
cmd_args+=("$par_output")
```

### 3. Conditional Parameter Inclusion

Use Bash parameter expansion for optional parameters:

```bash
# Include parameter only if variable is set and not empty
${par_threads:+--threads "$par_threads"}

# Include flag only if boolean is true (after unsetting false values)
${par_verbose:+--verbose}
```

### 4. Boolean Handling

Unset boolean parameters that are "false":

```bash
# Single parameter
[[ "$par_verbose" == "false" ]] && unset par_verbose

# Multiple parameters using loop
unset_if_false=(
    par_verbose
    par_quiet
    par_force
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done
```

### 5. Proper Quoting

Always quote variables that might contain spaces or special characters:

```bash
# Correct
--input "$par_input"
--output "$par_output"

# For special characters, use @Q expansion
--pattern "${par_pattern@Q}"
```

## Real-World Example

Here's an example from the bowtie2_build component:

```bash
#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_large_index" == "false" ]] && unset par_large_index
[[ "$par_noauto" == "false" ]] && unset par_noauto
[[ "$par_packed" == "false" ]] && unset par_packed

# Create output directory
mkdir -p "$par_output"

# Determine index basename
if [ -n "$par_index_name" ]; then
    index_basename="$par_index_name"
else
    index_basename=$(basename "$par_input" .fasta)
fi

# Build command arguments
cmd_args=(
    ${par_fasta:+-f}
    ${par_cmdline:+-c}
    ${par_large_index:+--large-index}
    ${par_noauto:+-a}
    ${par_packed:+-p}
    ${par_bmax:+--bmax "$par_bmax"}
    ${par_offrate:+-o "$par_offrate"}
    "$par_input"
    "$par_output/$index_basename"
)

# Execute bowtie2-build
bowtie2-build "${cmd_args[@]}"
```

## Advanced Patterns

### Multiple Input Handling

If your tool accepts multiple inputs with custom separators:

```bash
# Convert Viash's semicolon separator to comma
par_disable_filters=$(echo "$par_disable_filters" | tr ';' ',')

cmd_args=(
    --disable-filters "$par_disable_filters"
)
```

### Complex File Handling

```bash
# Ensure output directory exists
mkdir -p "$(dirname "$par_output")"

# Handle relative paths
input_path=$(realpath "$par_input")
output_path=$(realpath "$par_output")
```

### Resource Management

```bash
# Use available resources
cmd_args=(
    ${meta_cpus:+--threads "$meta_cpus"}
    ${meta_memory_mb:+--memory "${meta_memory_mb}M"}
)
```

## Common Pitfalls

### 1. Unquoted Variables
```bash
# Wrong - can break with spaces
cmd_args=(--input $par_input)

# Correct
cmd_args=(--input "$par_input")
```

### 2. Improper Boolean Handling
```bash
# Wrong - will include false booleans
cmd_args=(${par_verbose:+--verbose})

# Correct - unset false values first
[[ "$par_verbose" == "false" ]] && unset par_verbose
cmd_args=(${par_verbose:+--verbose})
```

### 3. Array Expansion
```bash
# Wrong - treats array as single string
tool $cmd_args

# Correct - expands array elements
tool "${cmd_args[@]}"
```

## Testing Your Script

Always test your script with:
- Empty/missing optional parameters
- Parameters with spaces
- Boolean true/false values
- Edge cases specific to your tool
