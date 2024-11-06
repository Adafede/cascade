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
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x) {
        p()
        if (nrow(x != 0)) {
          prepare_plot(dataframe = x)
        } else {
          NA
        }
      }
    )
}
