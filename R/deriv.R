#' Deriv
#'
#' @param x X
#' @param y Y
#'
#' @return The derivative
#'
#' @examples NULL
deriv <- function(x, y) {
  diff(y) / diff(x)
}
