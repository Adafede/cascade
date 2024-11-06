#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
.datatable.aware <- TRUE

.onLoad <- function(libname, pkgname) {
  strat <- ifelse(test = .Platform$OS.type == "unix",
    yes = "multicore",
    no = "sequential"
  )
  future::plan(strategy = strat)
  packageStartupMessage(
    "Running in ",
    strat,
    " with ",
    future::nbrOfWorkers(),
    " workers detected"
  )
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to  ", pkgname)
  message(format(utils::citation(pkgname)))
}
