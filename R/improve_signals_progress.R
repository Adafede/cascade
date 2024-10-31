#' Improve signals progress
#'
#' @include improve_signal.R
#'
#' @param xs XS
#'
#' @return A list of dataframes with improved signals
#'
#' @examples NULL
improve_signals_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      improve_signal(df = x |>
        dplyr::select(
          time,
          intensity
        ))
    }
  )
}
