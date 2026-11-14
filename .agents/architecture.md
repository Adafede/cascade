# Architecture & Data Contracts

## Module Flow (`DESCRIPTION` `Collate` Order)

1. **Metadata, data docs, palettes, globals**: `cascade-package.R` -> `data.R`
   -> `colors.R` -> `globals.R`
2. **I/O & caching (`load.R`)**:
   - `fread_gz()` caches decompressed `.csv.gz`/`.tsv.gz` tables under
     `tempdir()` using file path, size, and `mtime` in `cascade_env`; unlinks
     corrupt temp files on read failure.
   - `get_zenodo_file()` caches Zenodo snapshots under
     `tools::R_user_dir("cascade", "cache")`, keying LOTUS default under
     `"zenodo_lotus_path"` and custom `(doi, pattern)` downloads dynamically.
   - `load_annotations(show_example = TRUE)` returns a `data.table::copy()` with
     `mode` set to the requested polarity (`"pos"` or `"neg"`).
3. **Signal processing & peak fitting**:
   - `signal.R`: `runmean_center()` (ignores `NA`s), `filter_fft()`,
     `signal_sharpening()`, `improve_signal()` (finite grid `[rt_lo, rt_hi]`
     returning empty `data.table` when `< 5` points), `baseline_chromatogram()`,
     `normalize_chromatograms_list()`, `extract_chromatogram()` (resolves
     `.mzML` detectors by name first, then fallback index),
     `preprocess_chromatograms()`.
   - `get_peaks.R`: `find_peaks()` EGH/Gaussian solver + FWHM/area calculation.
   - `peaks.R`: `prepare_rt()`, `join_peaks()`, `normalize_chromato()`,
     `prepare_peaks()`, `extract_ms_progress()`, `compare_peaks()` (guarded via
     `safe_cor()` against zero-variance inputs).
4. **Annotations & plotting**:
   - `annotations.R`: `dt_left_join()`/`dt_inner_join()`/`dt_full_join()`,
     `keep_best_candidates()`, `prepare_comparison()` (`fixed = TRUE` literal
     organism matching), shared `encode_smiles()` helper.
   - `hierarchy.R`: base-R factor/rank helpers (`fct_reorder2_base()`,
     `dense_rank_base()` using `data.table` `.GRP`), `prepare_hierarchy()`,
     `prepare_plot()`, `prepare_plot_2()`, `treemaps_progress()`.
5. **Exported workflows**:
   - `check_chromatograms_alignment.R`
   - `check_peaks_integration.R`
   - `process_compare_peaks.R`
   - `generate_pseudochromatograms.R`
   - `generate_tables.R`
   - `query_wikidata.R` / `taxon_name_to_qid.R` / `generate_ids.R`
   - `prepare_tima_annotations.R` / `plot_tima.R`

## Critical Contracts --- Do Not Break

- **Detector argument names differ intentionally by function**:
  - `process_compare_peaks(type = ...)` uses `type` (`"baselined"`,
    `"improved"`, `"original"`).
  - `check_peaks_integration(chromatogram = ...)` uses `chromatogram`. Both call
    `match.arg(..., c("baselined", "improved", "original"))`.
- **Single-pass retention-time `shift` ownership**: Workflow callers
  (`check_peaks_integration()`, `process_compare_peaks()`,
  `generate_pseudochromatograms()`) pass `shift = 0` to
  `preprocess_chromatograms()`, because `preprocess_peaks()` -> `prepare_rt()`
  applies `shift` once across peaks. Calling `preprocess_chromatograms()`
  directly with `shift != 0` shifts all 6 returned tables (`df_xy` + `*_long`)
  uniformly.
- **Mode forwarding**: `generate_pseudochromatograms()` forwards `mode` (`"pos"`
  / `"neg"`) into `load_annotations()`; `generate_tables()` accepts
  `file_positive`, `file_negative`, or both (erroring if both are `NULL` when
  `show_example = FALSE`) and exports TSV (despite parameter `export_csv`).
- **`prepare_tima_annotations()` & `plot_tima()`**: unnamed annotation lists get
  integer names before appending `"_pred"`, accept `data.frame` entries, and
  `plot_tima()` merges multiple annotation tables with `rbindlist()`.
