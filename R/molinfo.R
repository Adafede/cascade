#' Molinfo
#'
#' @param x X
#'
#' @return A mol image
#'
#' @examples NULL
molinfo <- function(x) {
  gt::web_image(
    url = paste0(
      "https://www.simolecule.com/cdkdepict/depict/bot/svg?smi=",
      x
      # "&annotate=cip&r=0"
    ),
    height = as.numeric(75)
  )
}
