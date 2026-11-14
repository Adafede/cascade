# Testing Guide for `cascade`

`cascade` uses [`tinytest`](https://cran.r-project.org/package=tinytest) for
unit, regression, and offline integration testing. Test files live in
`inst/tinytest/`, and installed-package checks run via `tests/tinytest.R`.

## Test Suite Layout

  | File                                  | Scope                                                                                                                                                                                                                                                |
  | ---                                   | ---                                                                                                                                                                                                                                                  |
  | `inst/tinytest/test_cascade.R`        | End-to-end integration tests across the exported workflow functions using bundled examples and offline fixtures                                                                                                                                      |
  | `inst/tinytest/test_bugs.R`           | Regression tests for edge-case bugs in `compare_peaks()` and `get_peaks()`                                                                                                                                                                           |
  | `inst/tinytest/test_signal.R`         | Signal processing (`runmean_center`, `deriv`, `second_der`, `filter_fft`, `signal_sharpening`, `improve_signal`, `baseline_chromatogram`, `normalize_chromatograms_list`, `extract_chromatogram`, `preprocess_chromatograms`)                        |
  | `inst/tinytest/test_peaks.R`          | Peak integration, RT/m/z window preparation, MS1 XIC extraction, and `compare_peaks()`                                                                                                                                                               |
  | `inst/tinytest/test_get_peaks.R`      | Chromatographic peak detection (`find_peaks`, Gaussian/EGH/raw fitting, smoothing kernels, `sd.max`, time unit scaling)                                                                                                                              |
  | `inst/tinytest/test_load.R`           | Data loaders and cache management (`check_export_dir`, `fread_gz`, `get_zenodo_file`, `load_annotations`, `load_chromatograms`, `load_features*`, `load_ms_data`, `parse_ms_spectra`)                                                                |
  | `inst/tinytest/test_query_wikidata.R` | SPARQL client (`query_wikidata`), HTTP GET/POST, TSV quote/header handling, QLever failover, `queries_progress()`, `tables_progress()`, `format_gt()`, `molinfo()`, and `taxon_name_to_qid()`                                                        |
  | `inst/tinytest/test_annotations.R`    | Order-preserving `data.table` joins (`dt_left_join`, `dt_inner_join`, `dt_full_join`), `keep_best_candidates()`, `make_confident()`, `make_other()`, `no_other()`, `make_no_stereo()`, `make_chromatographiable()`, and `prepare_comparison()`       |
  | `inst/tinytest/test_hierarchy.R`      | Taxonomic/chemical hierarchy builders (`fct_reorder2_base`, `dense_rank_base`, `prepare_hierarchy`, `hierarchies_progress`, `hierarchies_grouped_progress`, `prepare_plot`, `prepare_plot_2`, `plot_results_2`, `treemaps_progress`)                 |
  | `inst/tinytest/test_workflows.R`      | High-level workflow branches (`check_chromatograms_alignment`, `check_peaks_integration`, `process_compare_peaks`, `generate_pseudochromatograms`, `generate_tables`, `generate_ids`, `prepare_tima_annotations`, `plot_tima`) with offline fixtures |

All tests run offline: network download tests isolate `R_USER_CACHE_DIR` in a
temporary directory and mock `curl::curl_download` with synthetic fixtures.

## Running Tests Locally

During development, load the working tree with `pkgload::load_all()` before
calling `tinytest::run_test_dir()` or `tinytest::run_test_file()`. Do **not**
use `tinytest::test_all(".")`, because `tinytest::test_all()` attaches an
installed copy of `cascade` if present in the R library rather than your
working-tree changes.

Run the full offline test suite against your working tree:

```r
pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
res <- tinytest::run_test_dir("inst/tinytest")
print(res)
stopifnot(all(as.logical(res)))
```

Run a single test file after loading the source tree:

```r
pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
res <- tinytest::run_test_file("inst/tinytest/test_signal.R")
print(res)
stopifnot(all(as.logical(res)))
```

Run tests against an already-installed package:

```r
tinytest::test_package("cascade")
```

Run a full tarball build and check (note that re-serializing ~45 MB of `bzip2`
lazy-load datasets under `data/` takes several minutes):

```bash
R CMD build .
R CMD check --no-manual cascade_*.tar.gz
```

## Static Checks and Spelling

```bash
air format . --check
jarl check .
panache format --check .
panache lint .
```

Run `lintr` and `spelling` in the same session after `pkgload::load_all()` so
that `object_usage_linter` inspects your current working tree rather than an
older installed library:

```r
pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
stopifnot(length(lintr::lint_package()) == 0L)
stopifnot(nrow(spelling::spell_check_package()) == 0L)
```

## Measuring Test Coverage

```r
cov <- covr::package_coverage(".")
print(cov)
covr::zero_coverage(cov)
```

## Testing Conventions

1. **Assert on concrete outcomes.** Verify exact dimensions, column names,
   values, types, and error messages rather than testing only that a call
   returns.
2. **Cover boundaries and failure paths.** Test empty inputs, single-row and
   short traces, `NA`/`NaN` values, zero-variance/flat signals, and invalid
   parameters via `expect_error(..., pattern = ...)` and
   `expect_warning(..., pattern = ...)`.
3. **Explicit numeric tolerances.** Use `expect_identical()` for integer,
   character, logical, and structural checks, and
   `expect_equal(..., tolerance = 1e-12)` for floating-point comparisons.
4. **State isolation.** `tinytest::run_test_file()` evaluates each top-level
   expression separately, so wrap temporary files, directories, env vars
   (`R_USER_CACHE_DIR`), and mocked bindings inside `local({ ... })` with
   `on.exit(..., add = TRUE)` or `withr::local_*()`. Never mutate the user's
   real package cache directory.
