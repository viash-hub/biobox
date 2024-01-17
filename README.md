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

- **Documentation of Functionality**: The purpose and functionality of
  each component should be adequately described.
- **Documentation of Inputs and Outputs**: All input and output
  arguments should have a description and example (with extension).
- **Docker Image**: A Docker image (with optional additional
  dependencies) should be provided.
- **Write unit tests**: A unit test with possibly test resources needs
  to be provided.
- **Provide test resources**: If the unit test requires test resources,
  these should be provided in the `test_resources` section of the
  component.
- **Versioning**: If the component uses custom software (not installed
  via Apt, Apk, Yum, Pip, Conda, or R), a Bash script `version.sh` needs
  to be provided that outputs the version of the software.
- **File format specifications**: If a component returns a directory or
  data structure such as AnnData or MuData, a specification of the file
  format should be provided.

See the [CONTRIBUTING](CONTRIBUTING.md) file for more details.

## Repository Structure

…

## Installation and Usage

…

## Support and Community

For support, questions, or to join our community:

- **Issues**: Submit questions or issues via the [GitHub issue
  tracker](https://github.com/openpipelines-bio/base/issues).
- **Discussions**: Join our discussions via [GitHub
  Discussions](https://github.com/openpipelines-bio/base/discussions).

## License

This repository is licensed under an MIT license. See the
[LICENSE](LICENSE) file for details.
