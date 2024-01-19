
# Contributing guidelines

We encourage contributions from the community. To contribute:

1. **Fork the Repository**: Start by forking this repository to your account.
2. **Develop Your Component**: Create your Viash component, ensuring it aligns with our best practices (detailed below).
3. **Submit a Pull Request**: After testing your component, submit a pull request for review.

## Documentation of Functionality

The purpose and functionality of each component should be adequately described.

Example:

```yaml
functionality:
  name: star_align
  namespace: bioinformatics
  description: |
    Aligns reads to a reference genome using STAR.
```

## Documentation of Inputs and Outputs

All input and output arguments should have a description and example (with extension).

```yaml
functionality:
  # ...
  arguments:
    - name: --input
      type: file
      description: Input reads in FASTQ format. If the file is compressed, it must have the extension `.gz`.
      example: input.fastq.gz
      required: true
    - name: --output
      type: file
      direction: output
      description: Output BAM file.
      example: output.bam
      required: true
```

## Docker Image

A Docker image (with optional additional dependencies) should be provided.

```yaml
functionality:
  # ...
platforms:
  - type: docker
    image: python:3.10
    setup:
      - type: python
        packages: numpy
```

This container should also have `ps` installed.

## Write unit tests

A unit test with possibly test resources needs to be provided.

```yaml
functionality:
  # ...
  test_resources:
    - type: python_script
      path: script.py
```

With `script.py`:

```python
# ... todo
```

The bare minimum of the unit test is to run the component and check whether the output exists. Ideally, the unit test should also check whether the output is correct.

## Provide test resources

If the unit test requires test resources, these should be provided in the `test_resources` section of the component.

```yaml
# ... todo
```

TODO: discuss hosting test resources

## Versioning

If the component uses custom software (not installed via Apt, Apk, Yum, Pip, Conda, or R), a Bash script `version.sh` needs to be provided that outputs the version of the software. 

The output of this script should be a yaml file with the version of each software as a string.

```yaml
functionality:
  # ...
  version:
    type: bash
    path: version.sh
```

With `version.sh`:

```bash
#!/bin/bash

cat <<-END_VERSIONS
star: "$(STAR --version | sed -e "s/STAR_//g")"
samtools: "$(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')"
gawk: "$(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*$//')"
END_VERSIONS
```

## File format specifications

If a component returns a directory or data structure such as AnnData or MuData, a specification of the file format should be provided.

### Directory file format specification

```yaml
functionality:
  # ...
  arguments:
    - name: --output
      type: file
      # ...
      example: output/
      info:
        format:
          - type: directory
            contents:
              - type: file
                name: counts.csv
                description: Normalised expression values
                required: true
              - type: file
                name: size_factors.csv
                description: The size factors created by the normalisation method, if any.
                required: false
```

### AnnData file format specification

```yaml
functionality:
  # ...
  arguments:
    - name: --output
      type: file
      # ...
      example: output.h5ad
      info:
        format:
          layers:
            - type: double
              name: normalized
              description: Normalised expression values
              required: true
          obs:
            - type: double
              name: size_factors
              description: The size factors created by the normalisation method, if any.
              required: false
          uns:
            - type: string
              name: normalization_id
              description: "Which normalization was used"
              required: true
```

## Workflow

### Step 1: Find a component to contribute

### Step 2: Add config template

Change all occurrences of `xxx` to the name of the component.

Contents of `src/xxx/config.vsh.yaml`:

```yaml
functionality:
  name: xxx
  description: xxx
  info:
    keywords: [tag1, tag2]
    homepage: yyy
    documentation: yyy
    repository: yyy
    reference: "doi:yyy"
    licence: yyy
  argument_groups:
    - name: Inputs
      arguments:
    - name: Outputs
      arguments:
    - name: Arguments
      arguments:
  resources:
    - type: bash_script
      path: script.sh
  test_resources:
    - type: bash_script
      path: test.sh
    - type: file
      path: test_data
platforms:
  - type: docker
    image: quay.io/biocontainers/xxx:0.1.0--py_0
    setup:
      - type: docker
        run: |
          echo "xxx: \"0.1.0\"" > /var/software_versions.txt
  - type: nextflow
```

### Step 3: Find container

Google `biocontainer xxx` and find the container that is most suitable. Typically the link will be `https://quay.io/repository/biocontainers/xxx?tab=tags`.

### Step 4: Create help file

```bash
docker run --rm -it -v `pwd`/src/xxx/:/xxx quay.io/biocontainers/xxx:tag
xxx --help > /xxx/help.txt
```

