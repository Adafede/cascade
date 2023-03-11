#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
transform_ms <- function(x) {
  feature <- seq_along(x)
  y <- mclapply(X = feature, function(z) {
    df <-
      data.frame(
        intensity = x[[z]][1]@intensity,
        rtime = x[[z]][1]@rtime
      ) |>
      dplyr::filter(!is.na(intensity)) |>
      dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
        min(intensity))) |>
      dplyr::filter(intensity >= 0.1) |>
      dplyr::mutate(rtime = (rtime - min(rtime)) / (max(rtime) -
        min(rtime))) |> ## see https://github.com/sneumann/xcms/issues/593
      dplyr::arrange(rtime)
    return(df)
  })
  return(y)
}
