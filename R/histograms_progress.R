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
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(x) {
        plot_histograms(
          dataframe = x,
          label = "Organism"
        )
      }
    )
}
