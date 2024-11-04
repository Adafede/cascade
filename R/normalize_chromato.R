#' Normalize chromato
#'
#' @param x X
#' @param list List
#'
#' @return A normalized chromato
#'
#' @examples NULL
normalize_chromato <- function(x, list = list) {
  list[[unique(x$id)]] |>
    tidytable::filter(time >= x$rt_min[1] &
      time <= x$rt_max[1]) |>
    tidytable::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
      min(intensity))) |>
    tidytable::filter(intensity >= 0.1)
}
