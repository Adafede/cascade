source(file = "R/middle_pts.R")
source(file = "R/deriv.R")

#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
second_der <- function(x, y) {
  deriv(middle_pts(x), deriv(x, y))
}
