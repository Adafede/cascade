#' Transform MS
#'
#' @param x X
#'
#' @return A list with transformed MS
#'
#' @examples NULL
transform_ms <- function(x) {
  custom_min <- function(x) {
    if (length(x) > 0) {
      min(x)
    } else {
      Inf
    }
  }
  custom_max <- function(x) {
    if (length(x) > 0) {
      max(x)
    } else {
      Inf
    }
  }
  purrr::map(.x = seq_along(x), .f = function(z, min_int = 0.1) {
    data.frame(
      intensity = x[z, 1]@intensity,
      rtime = x[z, 1]@rtime
    ) |>
      tidytable::filter(!is.na(intensity)) |>
      tidytable::mutate(intensity = (intensity - custom_min(intensity)) / (custom_max(intensity) -
        custom_min(intensity))) |>
      tidytable::filter(intensity >= min_int) |>
      ## see https://github.com/sneumann/xcms/issues/593
      tidytable::mutate(rtime = (rtime - custom_min(rtime)) / (custom_max(rtime) -
        custom_min(rtime))) |>
      tidytable::arrange(rtime)
  })
}
