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
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      plot_histograms(
        dataframe = x,
        label = "Organism"
      )
    }
  )
}
