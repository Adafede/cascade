#' Load features not informed
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#'
#' @return A table of non informed features
#'
#' @examples NULL
load_features_not_informed <- function(file = NULL, show_example = FALSE) {
  if (show_example) {
    utils::data(
      "cascade_features_not_informed",
      package = "cascade",
      envir = environment()
    )
    cascade_features_not_informed |>
      tidytable::tidytable()
  } else {
    file |>
      tidytable::fread()
  }
}
