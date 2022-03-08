#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
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
