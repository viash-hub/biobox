# Docker and Engine Best Practices

This guide covers best practices for setting up Docker engines and managing dependencies in biobox components.

## Table of Contents
- [Preferred Approach: Biocontainers](#preferred-approach-biocontainers)
- [Finding Biocontainers](#finding-biocontainers)
- [Version Detection](#version-detection)
- [Docker Run Syntax](#docker-run-syntax)
- [Custom Containers](#custom-containers)
- [Recommended Base Containers](#recommended-base-containers)
- [Multi-tool Containers](#multi-tool-containers)
- [Container Optimization](#container-optimization)
- [Testing Docker Setup](#testing-docker-setup)

## Preferred Approach: Biocontainers

### Basic Setup

```yaml
engines:
  - type: docker
    image: quay.io/biocontainers/bowtie2:2.5.4--he96a11b_6
    setup:
      - type: docker
        run:
          - bowtie2 --version 2>&1 | head -1 | sed 's/.*version /bowtie2: /' > /var/software_versions.txt
```

### Key Requirements

1. **Use specific versions**: Always pin to specific versions with build strings
2. **Include version detection**: Add setup commands to create `/var/software_versions.txt`
3. **Verify command availability**: Ensure the container has the required commands from `requirements.commands`

## Finding Biocontainers

### Search Strategy

1. **Google search**: `biocontainer <tool_name>`
2. **Direct URL**: `https://quay.io/repository/biocontainers/<tool_name>?tab=tags`
3. **Check version compatibility**: Choose the most recent stable version
4. **Verify build string**: Include the complete version tag with build string

### Version Selection

```yaml
# Good: Specific version with build string
image: quay.io/biocontainers/samtools:1.17--hd87286a_2

# Bad: Latest or incomplete version
image: quay.io/biocontainers/samtools:latest
image: quay.io/biocontainers/samtools:1.17
```

## Version Detection

### Common Patterns

```bash
# Pattern 1: Tool outputs "Tool version X.Y.Z"
tool --version 2>&1 | head -1 | sed 's/.*version /tool: /' > /var/software_versions.txt

# Pattern 2: Tool outputs just version number
echo "tool: $(tool --version 2>&1 | head -1)" > /var/software_versions.txt

# Pattern 3: Complex version output, extract numeric part
tool --version 2>&1 | grep -E "^[0-9]" | head -1 | sed 's/^/tool: /' > /var/software_versions.txt

# Pattern 4: Version in specific format
tool --version 2>&1 | awk '{print "tool: " $NF}' > /var/software_versions.txt
```

### Real Examples

```bash
# bowtie2
bowtie2 --version 2>&1 | head -1 | sed 's/.*version /bowtie2: /' > /var/software_versions.txt

# samtools
samtools --version 2>&1 | head -1 | sed 's/samtools /samtools: /' > /var/software_versions.txt

# fastqc
fastqc --version 2>&1 | sed 's/FastQC v/fastqc: /' > /var/software_versions.txt
```

### Testing Version Detection

Always test your version detection command:

```bash
# Test in the container
docker run quay.io/biocontainers/tool:version bash -c "
  tool --version 2>&1 | head -1 | sed 's/.*version /tool: /'
"
```

## Docker Run Syntax

### List vs Multiline Strings

**Preferred: List format**
```yaml
run:
  # Single commands
  - command1 arg1 arg2
  - command2 arg1 arg2
  # Chained commands
  - command1 && command2 && command3
```

**Alternative: Multiline strings (for complex commands)**
```yaml
run: |
  command1 arg1 arg2 && \
  command2 arg1 arg2 && \
  command3 arg1 arg2
```

**Important:** Comments inside multiline strings (`run: |`) become Dockerfile `RUN` commands and will break the build. Use comments before the `run:` key or use the list format.

## Custom Containers

### When to Use Custom Containers

Use custom containers when:
- No suitable biocontainer exists
- You need to install additional dependencies
- You need a specific base environment (R, Python, etc.)

### Python-based Tools

```yaml
engines:
  - type: docker
    image: python:3.10-slim
    setup:
      - type: apt
        packages: [wget, curl]
      - type: python
        packages: 
        - numpy~=x.x.x
        - pandas~=x.x.x
        - scipy~=x.x.x
      - type: docker
        run:
          - pip install your-tool
          - echo "your-tool: $(your-tool --version)" > /var/software_versions.txt
```

### R-based Tools

```yaml
engines:
  - type: docker
    image: rocker/r2u:24.04
    setup:
      - type: r
        cran: [devtools, BiocManager]
        bioc: [Biostrings, GenomicRanges]
      - type: docker
        run:
          - R --version | head -1 | sed 's/R version /R: /' > /var/software_versions.txt
```

### Compilation from Source

```yaml
engines:
  - type: docker
    image: ubuntu:22.04
    setup:
      - type: apt
        packages: [build-essential, cmake, git, wget]
      - type: docker
        run:
          - wget https://github.com/user/tool/archive/v1.0.tar.gz && tar -xzf v1.0.tar.gz
          - cd tool-1.0 && make && make install
          - echo "tool: 1.0" > /var/software_versions.txt
```

## Recommended Base Containers

### General Purpose
- **Ubuntu**: `ubuntu:22.04` - Good for compilation and apt packages
- **Alpine**: `alpine:latest` - Minimal size, apk packages
- **Debian**: `debian:bookworm-slim` - Stable, well-supported

### Language-Specific

#### Python
```yaml
# Basic Python
image: python:3.10-slim

# With scientific packages
image: python:3.10

# GPU-enabled
image: nvcr.io/nvidia/pytorch:23.08-py3
```

#### R
```yaml
# Fast package installation
image: rocker/r2u:24.04

# Tidyverse included
image: rocker/tidyverse:4.3.0

# Bioconductor base
image: bioconductor/bioconductor_docker:RELEASE_3_17
```

#### Node.js
```yaml
# LTS version
image: node:18-slim

# Alpine variant
image: node:18-alpine
```

#### Other Languages
```yaml
# Java
image: openjdk:11-jre-slim

# Go
image: golang:1.20-alpine

# Rust
image: rust:1.70-slim

# Ruby
image: ruby:3.1-slim
```

## Multi-tool Containers

### Installing Multiple Tools

```yaml
engines:
  - type: docker
    image: ubuntu:22.04
    setup:
      - type: apt
        packages: [wget, curl, build-essential]
      - type: docker
        run:
          # Install tool 1
          - wget https://tool1.com/download && install_tool1
          # Install tool 2  
          - wget https://tool2.com/download && install_tool2
          # Create version file
          - echo "tool1: $(tool1 --version)" > /var/software_versions.txt
          - echo "tool2: $(tool2 --version)" >> /var/software_versions.txt
```

## Container Optimization

### Layer Efficiency

```yaml
# Good: Combine related commands
setup:
  - type: docker
    run: |
      apt-get update && \
      apt-get install -y wget curl && \
      wget https://tool.com/download && \
      install_tool && \
      apt-get clean && \
      rm -rf /var/lib/apt/lists/*

# Bad: Separate layers for each command
setup:
  - type: apt
    packages: [wget, curl]
  - type: docker
    run: wget https://tool.com/download
  - type: docker
    run: install_tool
  - type: docker
    run: apt-get clean
```

## Testing Docker Setup

### Local Testing

```bash
# Build the container locally
viash build config.vsh.yaml --setup cachedbuild

# Test interactively
docker run -it biobox/namespace/component:latest bash

# Check installed tools
which tool
tool --version

# Verify version file
cat /var/software_versions.txt
```

### Viash Docker Debugging

```bash
# Inspect the generated Dockerfile
viash run config.vsh.yaml -- ---dockerfile

# Build with cached layers (faster)
viash run config.vsh.yaml -- ---setup cachedbuild ---verbose

# Build from scratch (clean build)
viash run config.vsh.yaml -- ---setup build ---verbose

# Enter interactive debugging session
viash run config.vsh.yaml -- ---debug
```

### Common Issues

1. **Command not found**: Tool not in PATH or not installed
2. **Version detection fails**: Command syntax varies between tools
3. **Permission issues**: Tools installed in wrong location
4. **Missing dependencies**: Tool requires additional libraries

### Debugging Commands

```bash
# Check what's installed
docker run container_name which tool
docker run container_name tool --help
docker run container_name ls -la /usr/local/bin/

# Check environment
docker run container_name env
docker run container_name echo $PATH
```
