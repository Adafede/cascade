#' Title
#'
#' @return
#' @export
#'
#' @examples
load_lotus <- function() {
  if (!file.exists(paths$inst$extdata$source$libraries$lotus)) {
    message("Downloading LOTUS")
    get_lotus(export = paths$inst$extdata$source$libraries$lotus)
  } else {
    message("LOTUS found")
  }
}
