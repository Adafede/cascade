require(package = gt, quietly = TRUE)

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
molinfo <- function(x) {
  gt::web_image(
    url = paste0("https://molinfo-de.nprod.net/molecule/smiles/", x, ".svg"),
    height = as.numeric(75)
  )
}
