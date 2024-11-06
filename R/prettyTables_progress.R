#' Pretty tables progress
#'
#' @param xs XS
#'
#' @return Pretty tables
#'
#' @examples NULL
prettyTables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x) {
        p()
        temp_gt_function(
          table = x,
          title = names(x),
          subtitle = "All compounds"
        )
      }
    )
}
