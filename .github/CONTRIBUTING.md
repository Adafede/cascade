# Contributing to cascade

This guide explains how to report issues and propose changes to `cascade`.

## Reporting bugs and proposing changes

Before starting work on a large change, please open an issue at
<https://github.com/adafede/cascade/issues> to describe the problem or design
goal. For bug reports, include a minimal reproducible example (using
`show_example = TRUE` when the issue reproduces with bundled data, or small
synthetic inputs when it depends on custom files) along with `sessionInfo()`.

## Setting up a development environment

1. Fork and clone the repository:

   ```bash
   git clone https://github.com/<your-username>/cascade.git
   cd cascade
   ```

2. Install package dependencies in R:

   ```r
   pak::local_install_dev_deps(".")
   ```

   On Linux and macOS, `mzR` requires NetCDF development libraries (`netcdf` via
   Homebrew on macOS; `libnetcdf-dev` on Debian/Ubuntu).

## Running checks locally

- Run the full offline test suite against your working tree without installing:

  ```r
  pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
  res <- tinytest::run_test_dir("inst/tinytest")
  stopifnot(all(as.logical(res)))
  ```

- Run formatting, lint, and spelling checks:

  ```bash
  air format . --check
  jarl check .
  panache format --check .
  panache lint .
  ```

  ```r
  pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
  stopifnot(length(lintr::lint_package()) == 0L)
  stopifnot(nrow(spelling::spell_check_package()) == 0L)
  ```

  See [TESTING.md](../TESTING.md) for single-file test execution and
  `R CMD check` instructions.

## Documentation and style

- Edit roxygen2 comments in `R/*.R` rather than editing `man/*.Rd` or
  `NAMESPACE` by hand, then regenerate help pages with `roxygen2::roxygenise()`.
- Format R code with `air format .` and keep lines within 80 characters.
- Document user-facing changes in `NEWS.md`.
- Write unit and regression tests in `inst/tinytest/` using `tinytest`. Keep
  tests offline and isolated inside `local({ ... })` with temporary
  cache/working directories.

## Code of Conduct

All contributors are expected to follow the [Code of
Conduct](CODE_OF_CONDUCT.md).
