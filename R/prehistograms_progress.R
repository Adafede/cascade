#' Prehistograms progress
#'
#' @include prepare_plot.R
#'
#' @param xs XS
#'
#' @return A list of prehistograms
#'
#' @examples NULL
prehistograms_progress <- function(xs) {
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(x) {
        if (nrow(x != 0)) {
          prepare_plot(dataframe = x)
        } else {
          NA
        }
      }
    )
}
