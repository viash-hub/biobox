
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
  - type: file
    path: test_data
engines:
  - <...>
runners:
  - type: executable
  - type: nextflow
```

### Step 3: Fill in the metadata

Fill in the relevant metadata fields in the config. Here is an example of the metadata of an existing component.

```yaml
name: arriba
description: Detect gene fusions from RNA-Seq data
keywords: [Gene fusion, RNA-Seq]
links:
  homepage: https://arriba.readthedocs.io/en/latest/
  documentation: https://arriba.readthedocs.io/en/latest/
  repository: https://github.com/suhrig/arriba
  issue_tracker: https://github.com/suhrig/arriba/issues
references:
  doi: 10.1101/gr.257246.119
  bibtex: |
    @article{
      ... a bibtex entry in case the doi is not available ...
    }
license: MIT
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


### Step 6: Create or fetch test data

To help develop the component, it's interesting to have some test data available. In most cases, we can use the test data from the Snakemake wrappers. 

To make sure we can reproduce the test data in the future, we store the command to fetch the test data in a file at `src/xxx/test_data/script.sh`.

```bash
cat <<EOF > src/xxx/test_data/script.sh

# clone repo
if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

# copy test data
cp -r /tmp/snakemake-wrappers/bio/xxx/test/* src/xxx/test_data
EOF
```

The test data should be suitable for testing this component. Ensure that the test data is small enough: ideally <1KB, preferably <10KB, if need be <100KB.

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

```yaml
engines:
  - type: docker
    image: quay.io/biocontainers/xxx:0.1.0--py_0
```

The container should have your tool installed, as well as `ps`.

If you didn't find a suitable container in the previous step, you can create a custom container. For example:

```yaml
engines:
  - type: docker
    image: python:3.10
    setup:
      - type: python
        packages: numpy
```

For more information on how to do this, see the [documentation](https://viash.io/guide/component/add-dependencies.html#steps-for-creating-a-custom-docker-platform).

Here is a list of base containers we can recommend:

* Bash: [`bash`](https://hub.docker.com/_/bash), [`ubuntu`](https://hub.docker.com/_/ubuntu)
* C#: [`ghcr.io/data-intuitive/dotnet-script`](https://github.com/data-intuitive/ghcr-dotnet-script/pkgs/container/dotnet-script)
* JavaScript: [`node`](https://hub.docker.com/_/node)
* Python: [`python`](https://hub.docker.com/_/python), [`nvcr.io/nvidia/pytorch`](https://catalog.ngc.nvidia.com/orgs/nvidia/containers/pytorch)
* R: [`eddelbuettel/r2u`](https://hub.docker.com/r/eddelbuettel/r2u), [`rocker/tidyverse`](https://hub.docker.com/r/rocker/tidyverse)
* Scala: [`sbtscala/scala-sbt`](https://hub.docker.com/r/sbtscala/scala-sbt)

### Step 11: Write a runner script

Next, we need to write a runner script that runs the tool with the input arguments. Create a Bash script named `src/xxx/script.sh` which runs the tool with the input arguments.

```bash
#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_option" == "false" ]] && unset par_option

xxx \
  --input "$par_input" \
  --output "$par_output" \
  ${par_option:+--option}
```

When building a Viash component, Viash will automatically replace the `## VIASH START` and `## VIASH END` lines (and anything in between) with environment variables based on the arguments specified in the config.

As an example, this is what the Bash script for the `arriba` component looks like:

```bash
#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_skip_duplicate_marking" == "false" ]] && unset par_skip_duplicate_marking
[[ "$par_extra_information" == "false" ]] && unset par_extra_information
[[ "$par_fill_gaps" == "false" ]] && unset par_fill_gaps

arriba \
  -x "$par_bam" \
  -a "$par_genome" \
  -g "$par_gene_annotation" \
  -o "$par_fusions" \
  ${par_known_fusions:+-k "${par_known_fusions}"} \
  ${par_blacklist:+-b "${par_blacklist}"} \
  # ...
  ${par_extra_information:+-X} \
  ${par_fill_gaps:+-I}
```

Notes:

* If your arguments can contain special variables (e.g. `$`), you can use quoting (need to find a documentation page for this) to make sure you can use the string as input. Example: `-x ${par_bam@Q}`.

* Optional arguments can be passed to the command conditionally using Bash [parameter expansion](https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html). For example: `${par_known_fusions:+-k ${par_known_fusions@Q}}`

* If your tool allows for multiple inputs using a separator other than `;` (which is the default Viash multiple separator), you can substitute these values with a command like: `par_disable_filters=$(echo $par_disable_filters | tr ';' ',')`.

* If you have a lot of boolean variables that you would like to unset when the value is `false`, you can avoid duplicate code by using the following syntax:

```bash
unset_if_false=(
    par_argument_1
    par_argument_2
    par_argument_3
    par_argument_4
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done
```

this code is equivalent to

```bash
[[ "$par_argument_1" == "false" ]] && unset par_argument_1
[[ "$par_argument_2" == "false" ]] && unset par_argument_2
[[ "$par_argument_3" == "false" ]] && unset par_argument_3
[[ "$par_argument_4" == "false" ]] && unset par_argument_4
```


### Step 12: Create test script

If the unit test requires test resources, these should be provided in the `test_resources` section of the component. 

```yaml
test_resources:
  - type: bash_script
    path: test.sh
  - type: file
    path: test_data
```

Create a test script at `src/xxx/test.sh` that runs the component with the test data. This script should run the component (available with `$meta_executable`) with the test data and check if the output is as expected. The script should exit with a non-zero exit code if the output is not as expected. For example:

```bash
#!/bin/bash

set -e

## VIASH START
## VIASH END

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_empty() {
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
assert_file_contains_regex() {
  grep -q -E "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains_regex() {
  grep -q -E "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
#############################################

echo "> Run $meta_name with test data"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/reads_R1.fastq" \
  --output "output.txt" \
  --option

echo ">> Check if output exists"
assert_file_exists "output.txt"

echo ">> Check if output is empty"
assert_file_not_empty "output.txt"

echo ">> Check if output is correct"
assert_file_contains "output.txt" "some expected output"

echo "> All tests succeeded!"
```

Notes:

* Do always check the contents of the output file. If the output is not deterministic, you can use regular expressions to check the output.

* If possible, generate your own test data instead of copying it from an external resource.

### Step 13: Create a `/var/software_versions.txt` file

For the sake of transparency and reproducibility, we require that the versions of the software used in the component are documented.

For now, this is managed by creating a file `/var/software_versions.txt` in the `setup` section of the Docker engine.

```yaml
engines:
  - type: docker
    image: quay.io/biocontainers/xxx:0.1.0--py_0
    setup:
      - type: docker
        # note: /var/software_versions.txt should contain:
        #   arriba: "2.4.0"
        run: |
          echo "xxx: \"0.1.0\"" > /var/software_versions.txt
```
