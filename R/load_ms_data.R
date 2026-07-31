#' Load MS data
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#'
#' @return MS data
#'
#' @examples NULL
load_ms_data <- function(file = NULL, show_example = FALSE) {
  if (show_example) {
    message(
      "Loading example MS data from package data, doing it on disk will be more efficient"
    )
    utils::data(
      list = "cascade_ms_data",
      package = "cascade",
      envir = environment()
    )
    cascade_ms_data
  } else {
    file |>
      MSnbase::readMSData(mode = "onDisk", msLevel. = 1)
  }
}
