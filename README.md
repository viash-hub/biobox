# Base repository for reusable Viash components


This repository is a collection of reproducible and reusable Viash
components.

## Objectives

- **Reusability**: Facilitating the use of components across various
  projects and contexts.
- **Reproducibility**: Guaranteeing that bioinformatics analyses can be
  reliably replicated.
- **Best Practices**: Adhering to established standards in software
  development and bioinformatics.

## Contributing

We encourage contributions from the community. To contribute:

1.  **Fork the Repository**: Start by forking this repository to your
    account.
2.  **Develop Your Component**: Create your Viash component, ensuring it
    aligns with our best practices (detailed below).
3.  **Submit a Pull Request**: After testing your component, submit a
    pull request for review.

## Contribution Guidelines

- **\### Step 1: Find a component to contribute**: \* Find a tool to
  contribute to this repo.
- **\### Step 2: Add config template**: Change all occurrences of `xxx`
  to the name of the component.
- **\### Step 3: Fill in the metadata**: Fill in the relevant metadata
  fields in the config. Here is an example of the metadata of an
  existing component.
- **\### Step 4: Find a suitable container**: Google
  `biocontainer <name of component>` and find the container that is most
  suitable. Typically the link will be
  `https://quay.io/repository/biocontainers/xxx?tab=tags`.
- **\### Step 5: Create help file**: To help develop the component, we
  store the `--help` output of the tool in a file at `src/xxx/help.txt`.
- **\### Step 6: Create or fetch test data**: To help develop the
  component, it’s interesting to have some test data available. In most
  cases, we can use the test data from the Snakemake wrappers.
- **\### Step 7: Add arguments for the input files**: By looking at the
  help file, we add the input arguments to the config file. Here is an
  example of the input arguments of an existing component.
- **\### Step 8: Add arguments for the output files**: By looking at the
  help file, we now also add output arguments to the config file.
- **\### Step 9: Add arguments for the other arguments**: Finally, add
  all other arguments to the config file. There are a few exceptions:
- **\### Step 10: Add a Docker engine**: To ensure reproducibility of
  components, we require that all components are run in a Docker
  container.
- **\### Step 11: Write a runner script**: Next, we need to write a
  runner script that runs the tool with the input arguments. Create a
  Bash script named `src/xxx/script.sh` which runs the tool with the
  input arguments.
- **\### Step 12: Create test script**:
- **\### Step 12: Create a `/var/software_versions.txt` file**: For the
  sake of transparency and reproducibility, we require that the versions
  of the software used in the component are documented.

See the [CONTRIBUTING](CONTRIBUTING.md) file for more details.

## Support and Community

For support, questions, or to join our community:

- **Issues**: Submit questions or issues via the [GitHub issue
  tracker](https://github.com/viash-hub/base/issues).
- **Discussions**: Join our discussions via [GitHub
  Discussions](https://github.com/viash-hub/base/discussions).

## License

This repository is licensed under an MIT license. See the
[LICENSE](LICENSE) file for details.
