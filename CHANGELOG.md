# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-10-12

### Added

- Protein multivariant reverse translation: HGVS allele input such as
  `p.[Arg175His;Gly245Asp]` via `--max-aa-changes`, producing the Cartesian
  product of per-position candidates as `c.`/`g.` allele notation.
- `--require-adjacent-aa-changes` to require multivariant positions to be
  consecutive.
- `CITATION.cff` for GitHub/Zenodo citation metadata, and a Citation section
  in the README.
- License section in the README.

### Fixed

- `--use-inv-notation` is now applied to `g.` output, not just `c.` output.

### Changed

- Copyright notice now names the specific author instead of "contributors".

## [1.1.0] - 2026-06-03

### Added

- Represent non-adjacent delins components as HGVS `g.` multi-variant
  alleles.
- Support transcript resolution from HGVS p. accessions.
- Allow HGVS-like input without a reference or `p.`/`c.`/`g.` prefix.

### Documentation

- Recommend an external tool for canonicalization and note limitations of
  `c.`-to-`g.` projection.

## [1.0.0] - 2026-03-18

### Added

- Batch execution function and README section for programmatic usage.
- `--skip` option for the reverse translation script.
- Optional output file for errors.
- Configurable CSV field length limit, for files with very wide fields.
- Option to normalize indel representation when comparing reverse
  translations, to avoid reporting differences in representation only.
- Instructions for using a local SeqRepo instance.

### Changed

- Switched from black to ruff for formatting and linting.
- `dotenv` is only loaded when invoked as a script, not when called
  programmatically.

### Documentation

- General project overview and documentation of supported protein variant
  types.

## [0.1.0] - 2026-03-17

Initial public release.

### Added

- Reverse-translation of protein HGVS variants to coding (`c.`) and genomic
  (`g.`) notation, including codon-local indel and inversion candidates.
- Option governing whether length-changing nucleotide alleles are allowed as
  reverse translations of stop-gain variants.
- Script for comparing reverse translations against reference data sets, with
  pass-through column support.
- License and packaging (`pyproject.toml`, `Makefile`).

[1.2.0]: https://github.com/bbi-lab/variant-translation/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/bbi-lab/variant-translation/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/bbi-lab/variant-translation/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/bbi-lab/variant-translation/releases/tag/v0.1.0
