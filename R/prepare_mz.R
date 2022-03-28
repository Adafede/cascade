#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
prepare_mz <- function(x) {
  feature <- seq_along(1:nrow(x))
  y <- mclapply(X = feature, function(z) {
    mzr <- x[z, ] |>
      dplyr::select(mzmin = mz_min, mzmax = mz_max) |>
      as.matrix()
    return(mzr)
  })
  return(y)
}
