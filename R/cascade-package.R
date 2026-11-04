#' @keywords internal
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0154 Small molecules
#' @concept EDAM:topic_3370 Analytical chemistry
#' @concept EDAM:topic_2258 Cheminformatics
#' @concept EDAM:topic_0780 Plant biology
#' @concept EDAM:operation_3803 Natural product identification
#' @concept EDAM:operation_3214 Spectral analysis
#' @concept EDAM:operation_3203 Chromatogram visualisation
#' @concept EDAM:operation_2942 Treemap visualisation
#' @concept EDAM:operation_3215 Peak detection
#' @concept EDAM:operation_3628 Chromatographic alignment
#' @concept EDAM:operation_3799 Quantification
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table := .I .N .SD as.data.table data.table fifelse fread
#' @importFrom data.table fwrite rbindlist set setcolorder setkey setnames
#' @importFrom data.table setorder uniqueN
## usethis namespace: end
NULL

.datatable.aware <- TRUE # nolint: object_name_linter.

.cascade_env <- new.env(parent = emptyenv())

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ", pkgname)
  packageStartupMessage(format(utils::citation(pkgname)))
}
