#' Load features
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#'
#' @return A table of features
#'
#' @examples NULL
load_features <- function(file = NULL,
                          show_example = FALSE) {
  if (show_example) {
    # feature_table |>
    #   saveRDS(file = "inst/extdata/features.rds")
    readRDS(system.file("extdata", "features.rds", package = "cascade"))
  } else {
    file |>
      tidytable::fread()
  }
}
