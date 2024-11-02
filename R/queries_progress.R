#' Queries progress
#'
#' @param xs XS
#' @param start Start
#' @param end End
#' @param limit Limit
#'
#' @return A list of queries
#'
#' @examples NULL
queries_progress <- function(xs,
                             start = "0",
                             end = "9999",
                             limit = "1000000") {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x, start, end, limit) {
      p()
      paste0(
        query_part_1,
        x,
        query_part_2,
        start,
        query_part_3,
        end,
        query_part_4,
        paste("\nLIMIT", limit)
      )
    },
    start = start,
    end = end,
    limit = limit
  )
}
