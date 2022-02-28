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
    test = !dir.exists(dirname(dirname(dir))),
    yes = dir.create(dirname(dirname(dir))),
    no = paste(dir, "exists")
  )
  ifelse(
    test = !dir.exists(dirname(dir)),
    yes = dir.create(dirname(dir)),
    no = paste(dir, "exists")
  )
  ifelse(
    test = !dir.exists(dir),
    yes = dir.create(dir),
    no = paste(dir, "exists")
  )
}
