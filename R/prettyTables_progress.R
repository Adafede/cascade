#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
prettyTables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      temp_gt_function(
        table = x,
        title = names(x),
        subtitle = "All compounds"
      )
    }
  )
}
