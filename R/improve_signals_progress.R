#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
improve_signals_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      improve_signal(df = x |>
        dplyr::select(time,
          intensity = UV.1_CAD_1_0
        ))
    }
  )
}
