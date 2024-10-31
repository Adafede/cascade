#' Second der
#'
#' @include middle_pts.R
#' @include deriv.R
#'
#' @param x X
#' @param y Y
#'
#' @return The second derivative
#'
#' @examples NULL
second_der <- function(x, y) {
  deriv(middle_pts(x), deriv(x, y))
}
