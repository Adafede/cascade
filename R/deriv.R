#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
deriv <- function(x, y) {
  diff(y) / diff(x)
}
