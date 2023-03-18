#' Title
#'
#' @return
#' @export
#'
#' @examples
load_lotus <- function() {
  if (!file.exists(paths$data$source$libraries$lotus)) {
    message("Downloading LOTUS")
    get_lotus(export = paths$data$source$libraries$lotus)
  } else {
    message("LOTUS found")
  }
}
