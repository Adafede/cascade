#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
.datatable.aware <- TRUE

.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ", pkgname)
  message(format(utils::citation(pkgname)))
}
