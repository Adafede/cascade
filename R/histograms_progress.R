#' Histograms progress
#'
#' @include plot_histograms.R
#'
#' @param xs XS
#'
#' @return A list of histograms
#'
#' @examples NULL
histograms_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x) {
        p()
        plot_histograms(
          dataframe = x,
          label = "Organism"
        )
      }
    )
}
