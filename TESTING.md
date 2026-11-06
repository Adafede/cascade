# Testing Guide for `cascade`

`cascade` uses [`tinytest`](https://cran.r-project.org/package=tinytest) for
unit, property-based, and integration testing. Test files live in
`inst/tinytest/` and are invoked via `tests/tinytest.R`.

## Test Suite Layout

  | File                                  | Scope                                                                                                                                                                                                                                          |
  | ------                                | -------                                                                                                                                                                                                                                        |
  | `inst/tinytest/test_cascade.R`        | End-to-end integration tests across the 5 exported workflow functions                                                                                                                                                                          |
  | `inst/tinytest/test_bugs.R`           | Regression tests for edge-case bugs in `compare_peaks()` and `get_peaks()`                                                                                                                                                                     |
  | `inst/tinytest/test_signal.R`         | Signal processing (`runmean_center`, `deriv`, `second_der`, `filter_fft`, `signal_sharpening`, `improve_signal`, `baseline_chromatogram`, `normalize_chromatograms_list`, `preprocess_chromatograms`)                                          |
  | `inst/tinytest/test_peaks.R`          | Peak integration, RT/m/z window preparation, MS1 XIC extraction, `closest_idx()`, and `compare_peaks()`                                                                                                                                        |
  | `inst/tinytest/test_get_peaks.R`      | Chromatographic peak detection (`find_peaks`, Gaussian/EGH/raw fitting, smoothing kernels, `sd.max`, time unit scaling)                                                                                                                        |
  | `inst/tinytest/test_load.R`           | Data loaders (`check_export_dir`, `fread_gz`, `get_zenodo_file`, `load_annotations`, `load_chromatograms`, `load_features*`, `load_ms_data`, `parse_ms_spectra`)                                                                               |
  | `inst/tinytest/test_query_wikidata.R` | SPARQL client (`query_wikidata`), HTTP GET/POST, CSV-to-TSV fallback, QLever failover, `queries_progress()`, `tables_progress()`, `format_gt()`, `molinfo()`, and `taxon_name_to_qid()`                                                        |
  | `inst/tinytest/test_annotations.R`    | Order-preserving `data.table` joins (`dt_left_join`, `dt_inner_join`, `dt_full_join`), `keep_best_candidates()`, `make_confident()`, `make_other()`, `no_other()`, `make_no_stereo()`, `make_chromatographiable()`, and `prepare_comparison()` |
  | `inst/tinytest/test_hierarchy.R`      | Taxonomic/chemical hierarchy builders (`fct_reorder2_base`, `dense_rank_base`, `prepare_hierarchy`, `hierarchies_progress`, `hierarchies_grouped_progress`, `prepare_plot`, `prepare_plot_2`, `plot_results_2`, `treemaps_progress`)           |
  | `inst/tinytest/test_workflows.R`      | High-level workflow branches (`check_chromatograms_alignment`, `prepare_tima_annotations`, `generate_tables`, `generate_ids`) with offline fixtures                                                                                            |

## Running Tests Locally

Run the full test suite:

```r
tinytest::test_all(".")
```

Run a single fast unit test file (<1s) during development:

```r
pkgload::load_all(".", quiet = TRUE)
tinytest::run_test_file("inst/tinytest/test_signal.R")
```

Run `R CMD check`:

```r
rcmdcheck::rcmdcheck(
  ".",
  args = c("--no-manual", "--ignore-vignettes"),
  build_args = c("--no-build-vignettes")
)
```

## Measuring Test Coverage

Compute line coverage per file and per function with `covr`:

```r
cov <- covr::package_coverage(".")
print(cov)
covr::zero_coverage(cov)
```

## Testing Conventions

1. **No coverage-only tests.** Every test must assert on concrete outcomes
   computed independently (exact values, dimensions, column names, types,
   attributes, and ordering), never by re-running the function under test.
2. **Boundary and equivalence-class coverage.** Test `0`, `1`, `n - 1`, `n`,
   `n + 1`, empty inputs (`numeric(0L)`, 0-row `data.table`s), `NA`, `NaN`,
   ties, duplicates, unsorted inputs, and exact threshold boundaries (`<` vs
   `<=`, `>` vs `>=`).
3. **Explicit numeric tolerances.** Use `expect_identical()` for integer,
   character, logical, and structural checks, and
   `expect_equal(..., tolerance = 1e-12)` for floating-point comparisons.
4. **Error and warning assertions.** Verify failure paths with
   `expect_error(..., pattern = ...)` and warnings with
   `expect_warning(..., pattern = ...)`.
5. **State isolation.** Scope temporary files, directories, and mocked bindings
   inside `local({ ... })` or `withr::with_dir(..., { ... })` with
   `on.exit(..., add = TRUE)`. Note that `tinytest::run_test_file()` evaluates
   each top-level expression separately, so `on.exit()` must always be placed
   inside a block or `local({})` rather than at top level.
6. **Bug-fix workflow.** When fixing a bug, commit a failing regression test
   first (`test:`), then commit the fix (`fix:`) separately.

## Running Mutation Checks

To verify that tests detect logic regressions (flipped comparison operators,
swapped boolean connectives, off-by-one index shifts, negated conditions, or
dropped branches), apply in-namespace mutations and run the corresponding
`tinytest` file:

```r
pkgload::load_all(".", quiet = TRUE)
ns <- asNamespace("cascade")
unlockBinding("runmean_center", ns)

# Example mutation: flip k <= 1L boundary to k < 1L
orig_fn <- get("runmean_center", envir = ns)
body_txt <- paste(deparse(body(orig_fn)), collapse = "\n")
mut_body <- parse(text = sub("k <= 1L", "k < 1L", body_txt, fixed = TRUE))
mut_fn <- orig_fn
body(mut_fn) <- mut_body[[1L]]
assign("runmean_center", mut_fn, envir = ns)

res <- tinytest::run_test_file("inst/tinytest/test_signal.R", verbose = 0)
stopifnot(!tinytest::all_pass(res)) # Mutant must be killed

# Always restore original binding
assign("runmean_center", orig_fn, envir = ns)
```
