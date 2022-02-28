require(package = purrr, quietly = TRUE)
require(package = yaml, quietly = TRUE)

source(file = "R/log_debug.R")

#' Title
#'
#' @return
#' @export
#'
#' @examples
parse_yaml_params <- function() {
  log_debug("Loading yaml parameters")
  suppressWarnings(params <-
    yaml::read_yaml(
      file = paths$config$default$file,
      handlers = list(
        seq = function(x) {
          purrr::flatten(x)
        }
      )
    ))
  if (file.exists(paths$config$params$file)) {
    suppressWarnings(params <-
      yaml::read_yaml(
        file = paths$config$params$file,
        handlers = list(
          seq = function(x) {
            purrr::flatten(x)
          }
        )
      ))
  }
  return(params)
}
