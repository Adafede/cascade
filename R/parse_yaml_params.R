#' Title
#'
#' @return Params
#'
#' @examples NULL
parse_yaml_params <- function() {
  message("Loading yaml parameters")
  suppressWarnings(params <-
    yaml::read_yaml(
      file = paths$params$default$file,
      handlers = list(
        seq = function(x) {
          purrr::flatten(x)
        }
      )
    ))
  if (file.exists(paths$params$user$file)) {
    suppressWarnings(params <-
      yaml::read_yaml(
        file = paths$params$user$file,
        handlers = list(
          seq = function(x) {
            purrr::flatten(x)
          }
        )
      ))
  }
  return(params)
}
