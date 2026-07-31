#' Load features
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#'
#' @return A table of features
#'
#' @examples NULL
load_features <- function(file = NULL, show_example = FALSE) {
  if (show_example) {
    utils::data("cascade_features", package = "cascade", envir = environment())
    cascade_features |>
      tidytable::tidytable()
  } else {
    file |>
      tidytable::fread()
  }
}
