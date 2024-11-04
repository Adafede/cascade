#' Load features not informed
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#'
#' @return A table of non informed features
#'
#' @examples NULL
load_features_not_informed <- function(file = NULL,
                                       show_example = FALSE) {
  if (show_example) {
    # features_not_informed |>
    #   saveRDS(file = "inst/extdata/features_not_informed.rds")
    readRDS(system.file("extdata", "features_not_informed.rds", package = "cascade"))
  } else {
    file |>
      tidytable::fread()
  }
}
