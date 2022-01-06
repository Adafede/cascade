#' Title
#'
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
check_export_dir <- function(dir) {
  ifelse(
    test = !dir.exists(dir),
    yes = dir.create(dir),
    no = paste(dir, "exists")
  )
}
