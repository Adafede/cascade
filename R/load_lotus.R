source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima-r/main/R/get_last_version_from_zenodo.R")

#' Title
#'
#' @return
#' @export
#'
#' @examples
load_lotus <- function() {
  if (!file.exists(paths$data$source$libraries$lotus)) {
    message("Downloading LOTUS")
    get_last_version_from_zenodo(
      doi = paths$url$lotus$doi,
      pattern = paths$urls$lotus$pattern,
      path = paths$data$source$libraries$lotus
    )
  } else {
    message("LOTUS found")
  }
}
