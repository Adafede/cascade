#' Transform MS
#'
#' @param x X
#'
#' @return A list with transformed MS
#'
#' @examples NULL
transform_ms <- function(x) {
  feature <- seq_along(x)
  future.apply::future_lapply(X = feature, function(z) {
    data.frame(
      intensity = x[[z]][1]@intensity,
      rtime = x[[z]][1]@rtime
    ) |>
      tidytable::filter(!is.na(intensity)) |>
      tidytable::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
        min(intensity))) |>
      tidytable::filter(intensity >= 0.1) |>
      ## see https://github.com/sneumann/xcms/issues/593
      tidytable::mutate(rtime = (rtime - min(rtime)) / (max(rtime) -
        min(rtime))) |>
      tidytable::arrange(rtime)
  })
}
