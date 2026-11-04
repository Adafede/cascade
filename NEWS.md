# cascade

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
