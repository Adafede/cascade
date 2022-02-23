#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
subtables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      x %>%
        dplyr::filter(chemical_pathway == .[1, "chemical_pathway"]) %>%
        dplyr::group_by(chemical_class) %>%
        dplyr::add_count(sort = TRUE) %>%
        dplyr::select(-n) %>%
        dplyr::group_by(chemical_superclass) %>%
        dplyr::add_count(sort = TRUE) %>%
        dplyr::select(-n, -chemical_pathway) %>%
        dplyr::distinct()
    }
  )
}
