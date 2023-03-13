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
  y <- mclapply(
    X = feature, FUN = function(z) {
      mzr <- x[z, ] |>
        select(mzmin = mz_min, mzmax = mz_max) |>
        as.matrix()
      return(mzr)
    }
  )
  return(y)
}
