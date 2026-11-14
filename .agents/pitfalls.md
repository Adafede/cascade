# Known Pitfalls

1. **Never use `tinytest::test_all(".")` during development.**
   `tinytest::test_all()` attaches an installed copy of `cascade` if one exists
   in the R library, silently ignoring edits in `R/` and failing or masking
   bugs. Always load the source tree first with:
   `pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)`
   then run `tinytest::run_test_dir("inst/tinytest")` or
   `tinytest::run_test_file(...)`.

2. **`lintr::lint_package()` requires `pkgload::load_all()` in the same
   session.** `object_usage_linter` inspects loaded namespaces. Without
   `load_all()` first, `lint_package()` checks against any old installed
   `cascade` and falsely flags new internal helpers (such as `encode_smiles`).

3. **`tinytest` evaluates top-level expressions individually.** Putting
   `on.exit()` or `withr::local_*()` at file top level executes cleanup after a
   single statement instead of at file end. Always wrap test bodies or blocks in
   `local({ ... })` and isolate `R_USER_CACHE_DIR` with a temporary directory
   (`~/Library/Caches/org.R-project.R/R/cascade` on macOS must remain untouched
   by tests).

4. **Anon closures in `assign()` trigger `lintr` inside `local()`.** When
   mocking namespace functions in `tinytest` inside `local({ ... })`, build
   closures with a factory helper (`force(val); function(...) ...`) rather than
   capturing outer `local()` bindings directly inside `assign()`.

5. **`.GRP` and `data.table` NSE globals.** `R/hierarchy.R` uses
   `data.table::.GRP` and other column symbols via NSE. Keep all used NSE
   symbols registered in `R/globals.R`; update `.lintr` or `R/globals.R`
   together. Always `data.table::copy()` bundled `cascade_*` objects before
   in-place `:=` mutation.

6. **Wikidata SPARQL TSV quirks.** Wikidata TSV exports quote fields containing
   commas/tabs (`"..."`) and escape inner quotes as `\"`. `query_wikidata()`
   passes `quote = ""` to `data.table::fread()` for TSV and unescapes `\"`, and
   strips leading `?` on 0-row TSV column headers. `queries_progress()`
   suppresses benign HTTP/socket `readLines()` warnings.

7. **`.onAttach()` during non-installed source loads.** Calling
   `utils::citation(package = pkgname, lib.loc = libname)` fails under
   `roxygen2` or when the working directory is not named `cascade` unless
   wrapped in `tryCatch(suppressWarnings(...), error = function(e) NULL)`.

8. **Pre-commit (`prek.toml`) rules.**
   - `.Rbuildignore` must remain `LC_ALL=C`-sorted (`file-contents-sorter`).
   - `.lintr`, `AGENTS.md`, `.agents/`, `.agent-work/`, `prek.toml` must remain
     in `.Rbuildignore` so `R CMD check` emits 0 hidden-file NOTEs.
   - `R/` must contain zero `print()`, `browser()`, or `debug()` statements.

9. **45 MB `bzip2` dataset latency.** `pkgload::load_all()` takes ~60 s and
   `R CMD INSTALL` / `R CMD check` takes ~10-15 min because `data/*.rda` are
   large `bzip2`-compressed datasets. Batch R checks into one session and avoid
   parallel heavy R processes.
