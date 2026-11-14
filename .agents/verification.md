# Verification Ladder & Definition of Done

## Verification Ladder (Cheapest First)

1. **Static & formatting (<2 s)**:

   ```bash
   air format . --check
   jarl check .
   panache format --check .
   panache lint .
   ```

2. **YAML & `.Rbuildignore` order (<1 s)**:

   ```bash
   LC_ALL=C sort -c .Rbuildignore
   ruby -ryaml -e 'Dir.glob(".github/**/*.y*ml").each { |f| YAML.safe_load_file(f) }'
   ```

3. **Targeted test file (~65 s)**:

   ```bash
   Rscript -e '
     pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
     r <- tinytest::run_test_file("inst/tinytest/<file>.R")
     print(r)
     stopifnot(all(as.logical(r)))
   '
   ```

4. **Full offline suite + `lintr` + `spelling` in one R session (~120 s)**:

   ```bash
   Rscript -e '
     pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
     l <- lintr::lint_package()
     print(l)
     stopifnot(length(l) == 0L)
     sp <- spelling::spell_check_package(vignettes = TRUE)
     print(sp)
     stopifnot(nrow(sp) == 0L)
     r <- tinytest::run_test_dir("inst/tinytest")
     print(r)
     stopifnot(all(as.logical(r)))
   '
   ```

5. **Pre-commit hooks**:

   ```bash
   prek run --all-files
   ```

6. **Clean-checkout tarball build & check (~10-15 min)**:

   ```bash
   git clone . /tmp/cascade_clean_verify
   cd /tmp/cascade_clean_verify
   R CMD build .
   _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual cascade_*.tar.gz
   ```

   Inspect tarball contents to verify `AGENTS.md`, `.agents/`, `.agent-work/`,
   `.lintr`, and `prek.toml` are absent, and check `00check.log` for
   `0 errors | 0 warnings | 0 notes`.

## Definition of Done

- All 490+ offline `tinytest` assertions pass with 0 unexpected warnings and no
  modification to user cache directories or working tree files.
- `air format . --check`, `jarl check .`, `panache format --check .`,
  `panache lint .`, `lintr::lint_package()`, and
  `spelling::spell_check_package(vignettes = TRUE)` all report 0 issues.
- If `@` roxygen comments changed, `roxygen2::roxygenise()` was run and updated
  `NAMESPACE` / `man/*.Rd` are staged in the same commit.
- Human-facing files remain free of emojis, filler adjectives, and references to
  agent files; `.agent-work/` stays gitignored and uncommitted.
