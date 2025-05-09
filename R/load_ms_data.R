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
    # ms_data <- file |>
    #   MSnbase::readMSData(mode = "inMemory", msLevel. = 1)
    message(
      "Loading example MS file in memory, doing it on disk will be more efficient"
    )
    # ms_data |>
    #   saveRDS(file = "inst/extdata/ms_data.rds")
    readRDS(system.file("extdata", "ms_data.rds", package = "cascade"))
  } else {
    file |>
      MSnbase::readMSData(mode = "onDisk", msLevel. = 1)
  }
}
