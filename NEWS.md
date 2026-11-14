# cascade

# cascade 0.0.0.9003

## Bug Fixes & Stability

- **Signal preprocessing and grid alignment (`R/signal.R`, `R/peaks.R`)**:
  - Updated `runmean_center()` to ignore `NA` values across rolling windows and
    handle interior `NA`s without zeroing surrounding signal.
  - Restricted `improve_signal()` interpolation grids to finite input bounds,
    returning an empty table when fewer than five points are present.
  - Resolved `.mzML` detectors by name in `extract_chromatogram()` before index
    fallback and added explicit validation in `check_peaks_integration()`,
    `process_compare_peaks()`, and `generate_pseudochromatograms()`.
  - Normalized retention-time `shift` handling so workflow callers apply `shift`
    once during peak preparation rather than double-shifting traces.
  - Added zero-variance and zero-range guards in `compare_peaks()`,
    `join_peaks()`, `normalize_chromato()`, and
    `normalize_chromatograms_list()`, and removed dead `closest_idx()`.
- **Workflows, tables, and caches (`R/annotations.R`, `R/generate_tables.R`,
  `R/generate_pseudochromatograms.R`, `R/load.R`,
  `R/prepare_tima_annotations.R`, `R/plot_tima.R`)**:
  - Forwarded `mode` to `load_annotations()` in `generate_pseudochromatograms()`
    and copied bundled example data before mutating columns.
  - Preserved negative-mode feature columns in `generate_tables()`, required at
    least one input file when `show_example = FALSE`, read LOTUS gzipped tables
    as character to avoid `melt()` coercion warnings, filtered `NA`/empty DOIs,
    and added literal matching (`fixed = TRUE`) in `prepare_comparison()`.
  - Fixed unnamed `annotations` labels (`"1_pred"` instead of `"NA_pred"`),
    `NULL` validation, `data.frame` inputs, and multi-sample binding in
    `prepare_tima_annotations()` and `plot_tima()`.
  - Included file size and modification time in `fread_gz()` cache keys, keyed
    non-LOTUS Zenodo downloads separately, and guarded `.onAttach()` citation
    lookup during source loading.
- **Wikidata SPARQL & example scripts (`R/query_wikidata.R`, `R/generate_ids.R`,
  `inst/scripts/`)**:
  - Fixed TSV quote parsing and `\"` unescaping, stripped leading `?` headers on
    0-row TSV responses, avoided duplicate SMILES URL-encoding via shared
    `encode_smiles()`, and added QID guards in `generate_ids()`.
  - Fixed misspelled filenames in `inst/scripts/`, removed commented dead code,
    directed outputs to `tempdir()`, and guarded live network calls in
    non-interactive runs.
- **Offline test isolation & static checks (`inst/tinytest/`, `.github/`)**:
  - Isolated all unit and workflow tests offline with temporary cache/working
    directories and mocked curl fixtures.
  - Added repository `.lintr`, enabled `spelling` errors and `codecov`
    enforcement, and repaired GitHub Actions workflows and package metadata.

## Mathematical & Algorithmic Improvements

- **Conjugate-symmetric Butterworth FFT low-pass filtering (`filter_fft()`)**:
  Replaced the asymmetric rectangular frequency-domain cutoff with a
  reflection-padded, conjugate-symmetric 4th-order Butterworth low-pass filter,
  reducing Gibbs phenomenon ringing and eliminating periodic boundary
  wrap-around artifacts.
- **Zero-phase derivative peak sharpening (`signal_sharpening()`,
  `runmean_symmetric()`, `improve_signal()`)**: Replaced the one-sided
  even-window moving average and offset finite-difference index slicing with
  zero-phase forward-backward smoothing and time-aligned centered 2nd and 4th
  derivatives on the interior grid `3:(n - 2)`, eliminating peak asymmetry
  distortion and retention time shift.
- **Automated cross-detector retention time shift estimation
  (`estimate_rt_shift()`, `check_chromatograms_alignment()`)**: Added sub-scan
  lag estimation via FFT cross-correlation (`cad_shift = "auto"` / `NULL` and
  `pda_shift = "auto"` / `NULL`).
- **Lan & Jorgenson EGH and Gaussian nonlinear least-squares peak fitting
  (`get_peaks()`)**: Implemented sub-grid half-height crossing interpolation,
  analytical Exponential-Gaussian Hybrid (EGH) parameter initialization
  ($\sigma$, $\tau$), bounded nonlinear least-squares refinement
  (`stats::nlminb()`), exact EGH $\text{FWHM}$ calculation, and trapezoidal area
  integration.
- **Trapezoidal peak integration and baseline clamping (`join_peaks()`,
  `trapz_integral()`, `preprocess_chromatograms()`)**: Replaced raw intensity
  summation with trapezoidal numerical quadrature over retention time and
  clamped negative baseline excursions to prevent global baseline inflation.
- **Window-preserving interpolated XIC-to-peak correlation
  (`normalize_chromato()`, `prepare_peaks()`, `transform_ms()`,
  `extract_ms_progress()`, `compare_peaks()`)**: Preserved physical peak
  retention-time window coordinates, recorded silent centroided MS1 scans as
  zero intensity rather than `NA`, and replaced nearest-neighbor snapping with
  continuous linear interpolation across the shared peak window.
- **Confidence-ordered 1D non-maximum suppression (`generate_tables()`)**:
  Replaced sequential lag-difference filtering (which dropped the
  earliest-eluting peak and suffered from chaining) with greedy 1D non-maximum
  suppression ordered by descending annotation confidence (`score_final`) and
  peak area within a `0.05` min retention-time radius.

# cascade 0.0.0.9002

## Architecture & Performance Overhaul

- Reduced external package dependencies from 17 to 7 (`baseline`, `curl`,
  `data.table`, `ggplot2`, `gt`, `mzR`, `plotly`), removing `BiocParallel`,
  `caTools`, `forcats`, `htmltools`, `httr2`, `MSnbase`, `purrr`, `R.utils`,
  `tidytable`, and `tima`.
- Consolidated 66 fragmented `R/*.R` files into 20 cohesive modules with
  explicit `Collate` ordering and full Roxygen2 documentation.
- Rewrote all tabular pipelines in idiomatic `data.table` (`fread()`,
  `fwrite()`, `foverlaps()`, `rbindlist()`, `melt()`, `dcast()`) and base R
  vectorization.
- Replaced `MSnbase::chromatogram()` with a vectorized `S4` MS1 scan extractor
  (`extract_ms_peak()`, `extract_ms_progress()`) that parses and binary-searches
  MS1 scans across all features simultaneously.
- Replaced `caTools::runmean()` with a fast base-R centered moving average
  (`runmean_center()`) and `MSnbase::closest()` with a vectorized binary-search
  matcher (`closest_idx()`).

## Completed In-Progress Functions & Palettes

- **Acetonitrile gradient response correction (`p_acn_i()`,
  `predict_response()`)**: Vectorized both functions in `R/annotations.R` and
  wired them into `prepare_comparison()` and `generate_pseudochromatograms()`
  via `correct_acn = FALSE`, `acn_time = 30.5`, and `dvol = 0.5`.
- **Signal-based sunburst plots (`plot_results_2()`)**: Completed in
  `R/hierarchy.R` using `data.table` and wired into
  `generate_pseudochromatograms()`, returning `plots_2`
  (`sunburst_signal_based_duo` and `sunburst_conf_signal_based_duo`) alongside
  `plots_1` and `treemaps`.
- **Interactive Wikidata/LOTUS `gt` tables (`format_gt()`)**: Completed in
  `R/query_wikidata.R` without requiring `htmltools` and wired into
  `generate_ids()` via `format_tables = FALSE` (returning `tables` alongside
  `plots`, `treemaps`, and `sunbursts`).
- **Titled treemaps and sunbursts (`treemaps_progress()`)**: Completed in
  `R/hierarchy.R` with explicit `tables` argument passing to display unique
  structure counts in plot titles, and wired into `generate_ids()` via
  `show_title = FALSE`.
- **Color palettes (`R/colors.R`)**: Preserved all qualitative and sequential
  color palettes (`microshades`, `microshades_grey`, `microshades_colors`,
  `paired`, `minibits`, `tableau20`, `discrete_rainbow_14`,
  `discrete_rainbow_23`, `elife`, `green_4`, `blue_4`, `pink_4`, `green_24`,
  `blue_24`, `pink_24`) and bundled them into `cascade_palettes`.

## Retired Prototype Functions

- `plot_histograms()`: Superseded by `plot_histograms_confident()` and
  `plot_histograms_taxo()` (for continuous retention-time pseudochromatograms)
  and `plot_histograms_litt()` (for discrete taxon/sample stacked bar charts).
- `plot_chromatogram()`: Unused 15-line `plotly` scratch wrapper expecting a
  non-standard `time` column (`cascade` uses `rtime`); superseded by
  `check_chromatograms_alignment()` and `plot_peak_detection()`.
- `prehistograms_progress()` and `histograms_progress()`: Trivial 1-line `map()`
  wrappers around `prepare_plot()` and `plot_histograms()`, replaced by direct
  `lapply()` calls in `generate_ids()` and `plot_tima()`.
- `save_histograms_progress()` and `save_treemaps_progress()`: Referenced
  external global pipeline variables (`paths`, `dimensions`, `size`) not present
  in the package and relied on deprecated `plotly::save_image()`.

## Ontology & Metadata

- Verified and standardized all EDAM ontology terms (`topic_3172`, `topic_0154`,
  `topic_3370`, `topic_2258`, `topic_0780`, `topic_0637`, `topic_0092`,
  `operation_3803`, `operation_3214`, `operation_3801`, `operation_3215`,
  `operation_3628`, `operation_3799`, `operation_3203`, `operation_2942`,
  `operation_3460`, `operation_3564`, `operation_0224`, `operation_0337`,
  `data_2536`, `data_0943`, `data_0962`, `data_1086`, `data_1463`, `data_1712`,
  `data_1872`, `format_3244`, `format_3475`, `format_3752`, `format_2331`,
  `format_3790`, `format_1196`, `format_1199`) against official `EDAM.tsv`
  across `DESCRIPTION`, `codemeta.json`, `CITATION.cff`, `inst/CITATION`, and
  `@concept` tags in `R/*.R`.

# cascade 0.0.9001

- Adapt to `tima` updates
- Added helper functions to visualize
  [TIMA](https://taxonomicallyinformedannotation.github.io/tima/) results
- Switched documentation from `pkgdown` to `altdoc`
- Updated minimal R version to `4.4.0` (and related Bioconductor dependencies)

# cascade 0.0.9000

- Initial version.
