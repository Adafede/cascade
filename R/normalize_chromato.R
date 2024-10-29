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
    dplyr::filter(time >= x$rt_min[1] &
      time <= x$rt_max[1]) |>
    dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
      min(intensity))) |>
    dplyr::filter(intensity >= 0.1)
}
