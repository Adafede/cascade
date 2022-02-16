#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
prehistograms_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      if (nrow(x != 0)) {
        prepare_plot(dataframe = x)
      } else {
        NA
      }
    }
  )
}
