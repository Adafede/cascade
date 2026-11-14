# Agent Guidelines for `cascade`

R package (`0.0.0.9003`, `R >= 4.5.0`, `LazyDataCompression: bzip2`) for
natural-extract LC-CAD/PDA-HRMS/MS analysis.

## Commands

Static checks (<2 s):
`air format . --check && jarl check . && panache format --check . && panache lint . && prek run --all-files`

Single `tinytest` file (~65 s due to ~45 MB `bzip2` dataset load):

```bash
Rscript -e 'pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE); r <- tinytest::run_test_file("inst/tinytest/test_signal.R"); print(r); stopifnot(all(as.logical(r)))'
```

Combined full offline tests + `lintr` + `spelling` in **one** R session (~120
s):

```bash
Rscript -e 'pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE); l <- lintr::lint_package(); print(l); stopifnot(length(l) == 0L); sp <- spelling::spell_check_package(vignettes = TRUE); print(sp); stopifnot(nrow(sp) == 0L); r <- tinytest::run_test_dir("inst/tinytest"); print(r); stopifnot(all(as.logical(r)))'
```

Roxygen regeneration (only when `@` tags change):
`Rscript -e 'roxygen2::roxygenise()'`

Full tarball build & check (~10-15 min):
`R CMD build . && R CMD check --no-manual cascade_*.tar.gz`

## Repo Map

- `R/`: 21 modules in `Collate` order (`data.table` + base R +
  `mzR`/`baseline`/`curl`/`ggplot2`/`gt`/`plotly`).
- `inst/tinytest/`: 10 offline test files (490 checks); isolate temp dirs and
  `R_USER_CACHE_DIR` with mocked `curl::curl_download`.
- `inst/scripts/`: 6 example scripts; live Wikidata/Zenodo scripts require
  `interactive()` or `CASCADE_RUN_LIVE=true`.
- `data/`: ~45 MB `bzip2` `.rda` example datasets. Do not alter schemas/values.
- `vignettes/` & `.github/workflows/`: Quarto articles (`cascade.qmd`, 4
  articles) and CI workflows.

## Invariants & Boundaries

- **Bundled data contract**: `cascade_features_informed` (`4132` rows) and
  `cascade_features_not_informed` (`650` rows) must remain reproduced by
  `process_compare_peaks(show_example = TRUE)`.
- **Generated files**: never hand-edit `NAMESPACE`, `man/*.Rd` (use
  `roxygen2::roxygenise()`), or `README.md` (edit `README.qmd` and run
  `quarto render README.qmd --to gfm`).
- **Two-surface separation**: human-facing files must never reference
  `AGENTS.md`, `CLAUDE.md`, or `.agents/` and must have zero emojis or marketing
  filler.
- Keep commits semantic (`fix`, `feat`, `refactor`, `test`, `docs`, `ci`,
  `chore`), subject `< 72` chars, no AI trailers, and never commit
  `.agent-work/`.

## Playbooks

- [`.agents/architecture.md`](.agents/architecture.md): module flow, parameter
  contracts, RT-shift ownership, caching.
- [`.agents/pitfalls.md`](.agents/pitfalls.md): `test_all` trap, `lintr`
  namespace trap, `tinytest` `on.exit`, Wikidata TSV quoting, `.Rbuildignore`
  sorting.
- [`.agents/verification.md`](.agents/verification.md): verification ladder and
  definition of done.
