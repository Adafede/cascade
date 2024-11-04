#' Subtables progress
#'
#' @param xs XS
#'
#' @return A list of subtables
#'
#' @examples NULL
subtables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      x |>
        dplyr::filter(chemical_pathway == .[1, "chemical_pathway"]) |>
        tidytable::group_by(chemical_class) |>
        tidytable::add_count(sort = TRUE) |>
        tidytable::select(-n) |>
        tidytable::group_by(chemical_superclass) |>
        tidytable::add_count(sort = TRUE) |>
        tidytable::select(-n, -chemical_pathway) |>
        tidytable::distinct()
    }
  )
}
