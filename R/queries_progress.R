#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
queries_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      paste0(
        query_part_1,
        x,
        query_part_2,
        params$references$date$start,
        query_part_3,
        params$references$date$end,
        query_part_4,
        paste("\nLIMIT", params$queries$limit)
      )
    }
  )
}
