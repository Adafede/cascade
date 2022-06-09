#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
baseline_chromatograms_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      baseline_chromatogram(df = x)
    }
  )
}