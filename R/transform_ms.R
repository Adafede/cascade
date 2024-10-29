#' Transform MS
#'
#' @param x X
#'
#' @return A list with transformed MS
#'
#' @examples NULL
transform_ms <- function(x) {
  feature <- seq_along(x)
  y <- future_lapply(X = feature, function(z) {
    df <-
      data.frame(
        intensity = x[[z]][1]@intensity,
        rtime = x[[z]][1]@rtime
      ) |>
      filter(!is.na(intensity)) |>
      mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
        min(intensity))) |>
      filter(intensity >= 0.1) |>
      mutate(rtime = (rtime - min(rtime)) / (max(rtime) -
        min(rtime))) |> ## see https://github.com/sneumann/xcms/issues/593
      arrange(rtime)
    return(df)
  })
  return(y)
}
