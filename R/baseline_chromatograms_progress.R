#' Baseline chromatograms progress
#'
#' @include baseline_chromatogram.R
#'
#' @param xs Xs
#'
#' @return A list of dataframes with baselined chromatogram
#'
#' @examples NULL
baseline_chromatograms_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      baseline_chromatogram(df = x)
    }
  )
}
