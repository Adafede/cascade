#' Title
#'
#' @param peak
#'
#' @return
#' @export
#'
#' @examples
extract_peak <- function(peak) {
  df <-
    data.frame(
      intensity = ms_chr[peak, 1]@intensity,
      time = ms_chr[peak, 1]@rtime
    ) |>
    dplyr::filter(!is.na(intensity))

  if (nrow(df) > 1) {
    f <- approxfun(
      x = df$time,
      y = df$intensity
    )

    timeow <- seq(
      from = min(df$time),
      to = max(df$time),
      by = (max(df$time) - min(df$time)) / 100
    )

    intensityeah <- f(seq(
      from = min(df$time),
      to = max(df$time),
      by = (max(df$time) - min(df$time)) / 100
    ))

    df_extrapolated <-
      data.frame(time = timeow, intensity = intensityeah)

    # plot(df_extrapolated)
    df_improved <-
      improve_signal(
        df = df_extrapolated,
        time_min = min(df$time),
        time_max = max(df$time)
      )
    # plot(df_improved)

    f <- approxfun(
      x = df_improved$time,
      y = df_improved$intensity
    )

    timeow <- seq(
      from = min(df$time),
      to = max(df$time),
      by = (max(df$time) - min(df$time)) / length(cad_ready$rtime)
    )

    intensityeah <- f(seq(
      from = min(df$time),
      to = max(df$time),
      by = (max(df$time) - min(df$time)) / length(cad_ready$rtime)
    ))

    ms_ready <-
      data.frame(rtime = timeow, intensity = intensityeah) |>
      dplyr::filter(!is.na(intensity)) |>
      dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
        min(intensity))) |>
      dplyr::filter(intensity >= 0.1) |>
      dplyr::mutate(rtime = (rtime - min(rtime)) / (max(rtime) -
        min(rtime))) |> #' see https://github.com/sneumann/xcms/issues/593
      dplyr::arrange(rtime)

    ms_peak <-
      MSnbase::Chromatogram(
        intensity = ms_ready$intensity,
        rtime = ms_ready$rtime
      )
  } else {
    ms_peak <- MSnbase::Chromatogram(intensity = 0, rtime = 0)
  }
  return(ms_peak)
}
