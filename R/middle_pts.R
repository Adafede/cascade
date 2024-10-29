#' Middle pts
#'
#' @param x X
#'
#' @return Middle pts
#'
#' @examples NULL
middle_pts <- function(x) {
  x[-1] - diff(x) / 2
}
