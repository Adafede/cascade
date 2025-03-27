#' Load features informed
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#'
#' @return A table of informed features
#'
#' @examples NULL
load_features_informed <- function(file = NULL, show_example = FALSE) {
  if (show_example) {
    # features_informed |>
    #   saveRDS(file = "inst/extdata/features_informed.rds")
    readRDS(system.file(
      "extdata",
      "features_informed.rds",
      package = "cascade"
    ))
  } else {
    file |>
      tidytable::fread()
  }
}
