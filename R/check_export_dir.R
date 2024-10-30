#' Check export dir
#'
#' @param dir Dir
#'
#' @return A log of checked dir
#'
#' @export
#'
#' @examples NULL
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
