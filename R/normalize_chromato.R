#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
normalize_chromato <- function(x) {
  chromatograms_cad_improved[[unique(x$id)]] |>
    dplyr::filter(time >= x$rt_min[1] &
      time <= x$rt_max[1]) |>
    dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
      min(intensity))) |>
    dplyr::filter(intensity >= 0.1)
}
