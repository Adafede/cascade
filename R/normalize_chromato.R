#' Normalize chromato
#'
#' @param x X
#' @param df_xy Df X Y
#'
#' @return A normalized chromato
#'
#' @examples NULL
normalize_chromato <- function(x, df_xy) {
  df_xy |>
    tidytable::filter(
      time >= x$rt_min[1] &
        time <= x$rt_max[1]
    ) |>
    tidytable::mutate(
      intensity = (intensity - min(intensity)) /
        (max(intensity) -
          min(intensity))
    ) |>
    tidytable::filter(intensity >= 0.1)
}
