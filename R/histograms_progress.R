#' Histograms progress
#'
#' @noRd
#'
#' @include plot_histograms.R
#'
#' @param xs XS
#'
#' @return A list of histograms
#'
#' @examples NULL
# jarl-ignore unused_function: <reason>
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
