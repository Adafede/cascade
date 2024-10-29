#' Prehistograms progress
#'
#' @param xs XS
#'
#' @return A list of prehistograms
#'
#' @examples NULL
prehistograms_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      if (nrow(x != 0)) {
        prepare_plot(dataframe = x)
      } else {
        NA
      }
    }
  )
}
