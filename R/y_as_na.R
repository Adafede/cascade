#' Y as NA
#'
#' @param x x
#' @param y y
#'
#' @return Y's replaced as NA's in X
#'
#' @examples NULL
y_as_na <- function(x, y) {
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  ## since ifelse wont work with factors
  ifelse(test = as.character(x) != y, yes = x, no = NA)
}
