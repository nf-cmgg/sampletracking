# nf-cmgg/sampletracking: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev

- Added the multiqc data to the output of the pipeline
- Fixed an issue where pool grouping for multiqc wasn't properly performed on pipeline resume
- Update modules
- Update to nf-core template 4.1.0
- Bump nf-schema to 3.0.0
- Revert the default configs base to the nf-core configs in preparation of converting nf-cmgg/configs to a private repo
- Add coverage based filter for crosscheck fingerprints
- Drop support for snp fastq inputs

## v1.0.2

- Fixed an issue where multiqc didn't run for each pool

## v1.0.1

- Bump modules
- Enable strict mode for Nextflow

## v1.0.0 - [21-01-2025]

Initial release of nf-cmgg/sampletracking, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
