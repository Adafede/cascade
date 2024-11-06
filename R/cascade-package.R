#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
.datatable.aware <- TRUE

.onLoad <- function(libname, pkgname) {
  ## Hack to avoid rcmdcheck warning since they are needed by {tidytable}
  ## to read gzip directly
  invisible(R.utils::capitalize("Bitter is better"))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ", pkgname)
  message(format(utils::citation(pkgname)))
}
