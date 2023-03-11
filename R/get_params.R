require(package = docopt, quietly = TRUE)
source(file = "r/parse_cli_params.R")
source(file = "r/parse_yaml_params.R")

#' Title
#'
#' @param step
#'
#' @return
#' @export
#'
#' @examples
get_params <- function(step) {
  doc_path <<- file.path(paths$inst$scripts$docopt, paste0(step, ".txt"))
  default_path <<-
    file.path(paths$params$default$path, paste0(step, ".yaml"))
  params_path <<-
    file.path(paths$params$user$path, paste0(step, ".yaml"))

  doc <<- readChar(
    con = doc_path,
    nchars = file.info(doc_path)$size
  )

  arguments <<- docopt(doc)

  params <<- parse_yaml_params()

  params <<- parse_cli_params()

  return(params)
}
