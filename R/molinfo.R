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
    url = paste0(
      "https://www.simolecule.com/cdkdepict/depict/bot/svg?smi=",
      x
      # "&annotate=cip&r=0"
    ),
    height = as.numeric(75)
  )
}
