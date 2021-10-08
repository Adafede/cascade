#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
middle_pts <- function(x) {
  x[-1] - diff(x) / 2
}
