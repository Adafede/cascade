#' Queries progress
#'
#' @param xs XS
#' @param start Start
#' @param end End
#' @param limit Limit
#' @param query_part_1 query_part_1
#' @param query_part_2 query_part_2
#' @param query_part_3 query_part_3
#' @param query_part_4 query_part_4
#'
#' @return A list of queries
#'
#' @examples NULL
queries_progress <- function(
  xs,
  start = "0",
  end = "9999",
  limit = "1000000",
  query_part_1,
  query_part_2,
  query_part_3,
  query_part_4
) {
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(
        x,
        start,
        end,
        limit,
        query_part_1,
        query_part_2,
        query_part_3,
        query_part_4
      ) {
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
      limit = limit,
      query_part_1 = query_part_1,
      query_part_2 = query_part_2,
      query_part_3 = query_part_3,
      query_part_4 = query_part_4
    )
}
