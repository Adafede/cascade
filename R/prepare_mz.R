#' Prepare mz
#'
#' @param x X
#'
#' @return A list of prepared mz's
#'
#' @examples NULL
prepare_mz <- function(x) {
  feature <- seq_along(1:nrow(x))
  y <- future.apply::future_lapply(
    X = feature, FUN = function(z) {
      mzr <- x[z, ] |>
        tidytable::select(mzmin = mz_min, mzmax = mz_max) |>
        as.matrix()
      return(mzr)
    }
  )
  return(y)
}
