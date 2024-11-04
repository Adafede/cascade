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
  packageStartupMessage("Welcome to CASCADE")
  strat <- ifelse(test = .Platform$OS.type == "unix",
    yes = "multicore",
    no = "sequential"
  )
  future::plan(strategy = strat)
  message("Running in ", strat, " with ", future::nbrOfWorkers(), " workers detected")
}
