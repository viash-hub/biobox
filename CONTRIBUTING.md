
# Contributing guidelines

We encourage contributions from the community. To contribute:

1. **Fork the Repository**: Start by forking this repository to your account.
2. **Develop Your Component**: Create your Viash component, ensuring it aligns with our best practices (detailed below).
3. **Submit a Pull Request**: After testing your component, submit a pull request for review.

## Procedure of adding a component

### Step 1: Find a component to contribute

* Find a tool to contribute to this repo.

* Check whether it is already in the [Project board](https://github.com/orgs/viash-hub/projects/1).

* Check whether there is a corresponding [Snakemake wrapper](https://github.com/snakemake/snakemake-wrappers/blob/master/bio) or [nf-core module](https://github.com/nf-core/modules/tree/master/modules/nf-core) which we can use as inspiration.

* Create an issue to show that you are working on this component.


### Step 2: Add config template

Change all occurrences of `xxx` to the name of the component.

Create a file at `src/xxx/config.vsh.yaml` with contents:

```yaml
name: xxx
description: xxx
keywords: [tag1, tag2]
links:
  homepage: yyy
  documentation: yyy
  issue_tracker: yyy
  repository: yyy
references: 
  doi: 12345/12345678.yz
license: MIT/Apache-2.0/GPL-3.0/...
requirements:
  commands: [command1, command2]
authors:
  - __merge__: /src/_authors/author_name.yaml
    roles: [author, maintainer]
argument_groups:
  - name: Inputs
    arguments: <...>
  - name: Outputs
    arguments: <...>
  - name: Arguments
    arguments: <...>
resources:
  - type: bash_script
    path: script.sh
test_resources:
  - type: bash_script
    path: test.sh
engines:
  - type: docker
    image: quay.io/biocontainers/xxx:version--build_string
    setup:
      - type: docker
        run: |
          xxx --version 2>&1 | head -1 | sed 's/.*version /xxx: /' > /var/software_versions.txt
runners:
  - type: executable
  - type: nextflow
```

### Step 3: Fill in the metadata

Fill in the relevant metadata fields in the config. Here is an example of the metadata of an existing component.

```yaml
name: bowtie2_build
namespace: bowtie2
description: |
  Build Bowtie2 index files from reference sequences.
  
  The Bowtie2 index is based on the FM Index of Ferragina and Manzini, which in turn 
  is based on the Burrows-Wheeler Transform. The algorithm used to build the index is 
  based on the blockwise algorithm of Karkkainen.
keywords: [Alignment, Indexing]
links:
  homepage: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  documentation: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
  repository: https://github.com/BenLangmead/bowtie2
references:
  doi: 10.1038/nmeth.1923
license: GPL-3.0
requirements:
  commands: [bowtie2-build]
authors:
  - __merge__: /src/_authors/robrecht_cannoodt.yaml
    roles: [author, maintainer]
```

**Important notes:**

* **Use `requirements.commands`**: Always specify the commands that your component requires in the `requirements.commands` field. This helps document dependencies and enables better validation.

* **Use `__merge__` for authors**: Instead of specifying author details inline, use the `__merge__` syntax to reference author information from the `/src/_authors/` directory. This promotes consistency and maintainability.

* **Author roles**: Specify appropriate roles for each author (`author`, `maintainer`, `contributor`, etc.).

### Step 3.1: Specify requirements

The `requirements` section documents the dependencies needed by your component:

```yaml
requirements:
  commands: [bowtie2-build, bowtie2]
```

**Why specify commands:**
- Documents which executables the component expects
- Enables validation that the Docker container has required tools
- Helps users understand dependencies
- Facilitates automated testing and CI/CD

**Examples:**
```yaml
# Single command
requirements:
  commands: [samtools]

# Multiple commands
requirements:
  commands: [bwa, samtools, bgzip]

# Commands with different names than the component
requirements:
  commands: [bowtie2-build]  # for bowtie2_build component
```

### Step 4: Find a suitable container

Google `biocontainer <name of component>` and find the container that is most suitable. Typically the link will be `https://quay.io/repository/biocontainers/xxx?tab=tags`.

If no such container is found, you can create a custom container in the next step. 


### Step 5: Create help file

To help develop the component, we store the `--help` output of the tool in a file at `src/xxx/help.txt`.

````bash
cat <<EOF > src/xxx/help.txt
```sh
xxx --help
```
EOF

docker run quay.io/biocontainers/xxx:tag xxx --help >> src/xxx/help.txt
````

Notes:

* This help file has no functional purpose, but it is useful for the developer to see the help output of the tool.

* Some tools might not have a `--help` argument but instead have a `-h` argument. For example, for `arriba`, the help message is obtained by running `arriba -h`:
  
  ```bash
  docker run quay.io/biocontainers/arriba:2.4.0--h0033a41_2 arriba -h
  ```


### Step 6: Create test data

**Best Practice: Generate test data in the test script rather than storing it in the repository.**

Instead of storing test data files in the repository or fetching them from external sources, it's preferred to generate test data programmatically within the test script. This approach:

- Keeps the repository size small
- Ensures reproducibility
- Makes tests self-contained
- Avoids external dependencies

For example, create a function in your test script to generate minimal test data:

```bash
# --- Helper function to create test FASTA ---
create_test_fasta() {
  file_path="$1"
  
  cat << 'EOF' > "$file_path"
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
EOF
}

# Usage in test
create_test_fasta "$TEMP_DIR/test_ref.fasta"
```

**When to use external test data:**

Only use external test data if:
- The tool requires very specific file formats that are difficult to generate
- You need real-world data to validate complex algorithms
- The test data is very small (<1KB)

If you must use external test data, prefer copying from established sources like Snakemake wrappers:

```bash
# If absolutely necessary, document how to fetch test data
cat <<EOF > src/xxx/test_data/script.sh
# clone repo
if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

# copy test data
cp -r /tmp/snakemake-wrappers/bio/xxx/test/* src/xxx/test_data
EOF
```

The test data should be suitable for testing this component. Ensure that external test data is small enough: ideally <1KB, preferably <10KB, if need be <100KB.

### Step 7: Add arguments for the input files

By looking at the help file, we add the input arguments to the config file. Here is an example of the input arguments of an existing component.

For instance, in the [arriba help file](src/arriba/help.txt), we see the following:

    Usage: arriba [-c Chimeric.out.sam] -x Aligned.out.bam \
                  -g annotation.gtf -a assembly.fa [-b blacklists.tsv] [-k known_fusions.tsv] \
                  [-t tags.tsv] [-p protein_domains.gff3] [-d structural_variants_from_WGS.tsv] \
                  -o fusions.tsv [-O fusions.discarded.tsv] \
                  [OPTIONS]

    -x FILE  File in SAM/BAM/CRAM format with main alignments as generated by STAR 
              (Aligned.out.sam). Arriba extracts candidate reads from this file. 

Based on this information, we can add the following input arguments to the config file.

```yaml
argument_groups:
  - name: Inputs
    arguments:
    - name: --bam
      alternatives: -x
      type: file
      description: |
        File in SAM/BAM/CRAM format with main alignments as generated by STAR
        (`Aligned.out.sam`). Arriba extracts candidate reads from this file.
      required: true
      example: Aligned.out.bam
```

Check the [documentation](https://viash.io/reference/config/functionality/arguments) for more information on the format of input arguments.

Several notes:

* Argument names should be formatted in `--snake_case`. This means arguments like `--foo-bar` should be formatted as `--foo_bar`, and short arguments like `-f` should receive a longer name like `--foo`.

* Input arguments can have `multiple: true` to allow the user to specify multiple files.

* The description should be formatted in markdown.

### Step 8: Add arguments for the output files

By looking at the help file, we now also add output arguments to the config file.

For example, in the [arriba help file](src/arriba/help.txt), we see the following:


    Usage: arriba [-c Chimeric.out.sam] -x Aligned.out.bam \
                  -g annotation.gtf -a assembly.fa [-b blacklists.tsv] [-k known_fusions.tsv] \
                  [-t tags.tsv] [-p protein_domains.gff3] [-d structural_variants_from_WGS.tsv] \
                  -o fusions.tsv [-O fusions.discarded.tsv] \
                  [OPTIONS]

     -o FILE  Output file with fusions that have passed all filters. 

     -O FILE  Output file with fusions that were discarded due to filtering. 

Based on this information, we can add the following output arguments to the config file.

```yaml
argument_groups:
  - name: Outputs
    arguments:
      - name: --fusions
        alternatives: -o
        type: file
        direction: output
        description: |
          Output file with fusions that have passed all filters.
        required: true
        example: fusions.tsv
      - name: --fusions_discarded
        alternatives: -O
        type: file
        direction: output
        description: |
          Output file with fusions that were discarded due to filtering. 
        required: false
        example: fusions.discarded.tsv
```

Note: 

* Preferably, these outputs should not be directories but files. For example, if a tool outputs a directory `foo/` containing files `foo/bar.txt` and `foo/baz.txt`, there should be two output arguments `--bar` and `--baz` (as opposed to one output argument which outputs the whole `foo/` directory).

### Step 9: Add arguments for the other arguments

Finally, add all other arguments to the config file. There are a few exceptions:

* Arguments related to specifying CPU and memory requirements are handled separately and should not be added to the config file.

* Arguments related to printing the information such as printing the version (`-v`, `--version`) or printing the help (`-h`, `--help`) should not be added to the config file.

* If the help lists defaults, do not add them as defaults but to the description. Example: `description: <Explanation of parameter>. Default: 10.`

Note:
  
* Prefer using `boolean_true` over `boolean_false`. This avoids confusion when specifying values for this argument in a Nextflow workflow.
  For example, consider the CLI option `--no-indels` for `cutadapt`. If the config for `cutadapt` would specify an argument `no_indels` of type `boolean_false`,
  the script of the component must pass a `--no-indels` argument to `cutadapt` when `par_no_indels` is set to `false`. This becomes problematic setting a value for this argument using `fromState` in a nextflow workflow: with `fromState: ["no_indels": true]`, the value that gets passed to the script is `true` and the `--no-indels` flag would *not* be added to the options for `cutadapt`. This is inconsitent to what one might expect when interpreting `["no_indels": true]`.
  When using `boolean_true`, the reasoning becomes simpler because its value no longer represents the effect of the argument, but wether or not the flag is set.

### Step 10: Add a Docker engine

To ensure reproducibility of components, we require that all components are run in a Docker container.

**Preferred approach - Use biocontainers:**

```yaml
engines:
  - type: docker
    image: quay.io/biocontainers/xxx:2.5.4--he96a11b_6
    setup:
      - type: docker
        run: |
          xxx --version 2>&1 | head -1 | sed 's/.*version /xxx: /' > /var/software_versions.txt
```

**Key requirements:**

1. **Use specific versions**: Always pin to specific versions with build strings (e.g., `2.5.4--he96a11b_6`)
2. **Include version detection**: Add setup commands to create `/var/software_versions.txt`
3. **Verify command availability**: Ensure the container has the required commands from `requirements.commands`

**If no biocontainer exists, create a custom container:**

```yaml
engines:
  - type: docker
    image: python:3.10
    setup:
      - type: python
        packages: [numpy, scipy]
      - type: docker
        run: |
          python --version | sed 's/Python /python: /' > /var/software_versions.txt
```

**Recommended base containers:**

* **Bash**: [`bash`](https://hub.docker.com/_/bash), [`ubuntu`](https://hub.docker.com/_/ubuntu)
* **C#**: [`ghcr.io/data-intuitive/dotnet-script`](https://github.com/data-intuitive/ghcr-dotnet-script/pkgs/container/dotnet-script)
* **JavaScript**: [`node`](https://hub.docker.com/_/node)
* **Python**: [`python`](https://hub.docker.com/_/python), [`nvcr.io/nvidia/pytorch`](https://catalog.ngc.nvidia.com/orgs/nvidia/containers/pytorch)
* **R**: [`eddelbuettel/r2u`](https://hub.docker.com/r/eddelbuettel/r2u), [`rocker/tidyverse`](https://hub.docker.com/r/rocker/tidyverse)
* **Scala**: [`sbtscala/scala-sbt`](https://hub.docker.com/r/sbtscala/scala-sbt)

For more information on custom containers, see the [Viash documentation](https://viash.io/guide/component/add-dependencies.html#steps-for-creating-a-custom-docker-platform).

### Step 11: Write a runner script

Next, we need to write a runner script that runs the tool with the input arguments. Create a Bash script named `src/xxx/script.sh` which runs the tool with the input arguments.

**Modern script structure using array-based arguments:**

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

**Key improvements in modern scripts:**

1. **Use `set -eo pipefail`**: Ensures the script fails fast on errors
2. **Array-based arguments**: Use a single `cmd_args` array instead of repetitive `cmd_args+=()` calls
3. **Conditional parameter inclusion**: Use `${var:+--flag "$var"}` for optional parameters
4. **Proper quoting**: Use `"$var"` for variables that might contain spaces

**Example from a real component:**

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

**Notes:**

* If your arguments can contain special variables (e.g. `$`), you can use `${par_bam@Q}` for proper escaping.

* Optional arguments can be passed conditionally using Bash [parameter expansion](https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html): `${par_known_fusions:+-k "$par_known_fusions"}`

* If your tool allows multiple inputs with custom separators, use: `par_disable_filters=$(echo $par_disable_filters | tr ';' ',')`

* For many boolean variables, you can use a loop to avoid repetition:

```bash
unset_if_false=(
    par_argument_1
    par_argument_2
    par_argument_3
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done
```


### Step 12: Create test script

**Best Practice: Generate test data within the test script instead of including static test resources.**

The test script should be self-contained and generate its own test data rather than relying on external files. This keeps the repository lean and ensures reproducibility.

```yaml
test_resources:
  - type: bash_script
    path: test.sh
```

Create a test script at `src/xxx/test.sh` that:
1. Generates its own test data
2. Runs the component with the generated data
3. Validates the output
4. Uses proper error handling

**Example test script structure:**

```bash
#!/bin/bash

set -e

TEMP_DIR="$meta_temp_dir"

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
create_test_fasta "$TEMP_DIR/input.fasta"

echo ">> Running $meta_name..."
"$meta_executable" \
  --input "$TEMP_DIR/input.fasta" \
  --output "$TEMP_DIR/output"

echo ">> Checking output exists..."
assert_dir_exists "$TEMP_DIR/output"

echo ">> Checking output is not empty..."
assert_file_not_empty "$TEMP_DIR/output/some_file"

echo ">> Checking output content..."
assert_file_contains "$TEMP_DIR/output/some_file" "expected_pattern"

# --- Test Case 2: Edge cases ---
echo ">>> Test 2: Edge case testing"
# ... additional test cases

echo "> All tests succeeded!"
```

**Key principles for test scripts:**

1. **Generate test data**: Use functions to create minimal test data rather than storing files
2. **Use `$meta_temp_dir`**: All temporary files should go in the provided temp directory
3. **Test multiple scenarios**: Include basic functionality, edge cases, and error conditions
4. **Validate outputs**: Always check that outputs exist, have expected content, and are not empty
5. **Use helper functions**: Include assertion functions for common checks
6. **Clear test structure**: Use descriptive echo statements to show test progress
7. **Fail fast**: Use `set -e` and proper assertions that exit on failure

**When you might need static test resources:**

Only include static test files in `test_resources` if:
- The tool requires very specific, complex file formats
- Generating equivalent test data is impractical
- You need real-world data to validate complex algorithms

If you must use static test resources:

```yaml
test_resources:
  - type: bash_script
    path: test.sh
  - type: file
    path: test_data
```

### Step 13: Create a `/var/software_versions.txt` file

For transparency and reproducibility, we require that software versions are documented in a standardized format.

This is managed by adding version detection commands to the `setup` section of the Docker engine:

```yaml
engines:
  - type: docker
    image: quay.io/biocontainers/xxx:2.5.4--he96a11b_6
    setup:
      - type: docker
        run: |
          xxx --version 2>&1 | head -1 | sed 's/.*version /xxx: /' > /var/software_versions.txt
```

**Best practices for version detection:**

1. **Use the most reliable version command**: Many tools support `--version`, some use `-v` or `version`
2. **Handle stderr output**: Some tools output version info to stderr, so use `2>&1`
3. **Extract clean version strings**: Use `sed` to create clean "tool: version" format
4. **Test the version command**: Verify that your version extraction works with the specific tool

**Common version extraction patterns:**

```bash
# For tools that output "Tool version X.Y.Z"
tool --version 2>&1 | head -1 | sed 's/.*version /tool: /' > /var/software_versions.txt

# For tools that output just the version number
echo "tool: $(tool --version 2>&1 | head -1)" > /var/software_versions.txt

# For tools with complex version output
tool --version 2>&1 | grep -E "^[0-9]" | head -1 | sed 's/^/tool: /' > /var/software_versions.txt
```

**Example from bowtie2:**

```yaml
setup:
  - type: docker
    run: |
      bowtie2-build --version 2>&1 | head -1 | sed 's/.*version /bowtie2-build: /' > /var/software_versions.txt
```

This ensures that the final `/var/software_versions.txt` contains entries like:
```
bowtie2-build: 2.5.4
```
