#' Prepare mz
#'
#' @param x X
#'
#' @return A list of prepared mz's
#'
#' @examples NULL
prepare_mz <- function(x) {
  feature <- seq_along(seq_len(nrow(x)))
  future.apply::future_lapply(
    X = feature,
    FUN = function(z) {
      x[z, ] |>
        tidytable::select(mzmin = mz_min, mzmax = mz_max) |>
        as.matrix()
    }
  )
}
