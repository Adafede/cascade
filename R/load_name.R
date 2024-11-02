#' Load name
#'
#' @param file File
#' @param default Default
#' @param show_example Show example? Default to FALSE
#'
#' @return A name
#'
#' @examples NULL
load_name <- function(file = NULL,
                      default = "210619_AR_06_V_03_2_01.mzML",
                      show_example = FALSE) {
  if (!show_example) {
    file |>
      basename()
  } else {
    default
  }
}
