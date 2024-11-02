#' Improve signal
#'
#' @include filter_fft.R
#' @include signal_sharpening.R
#'
#' @param df Dataframe
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param resample Resample
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A dataframe with improved signal
#'
#' @examples NULL
improve_signal <-
  function(df,
           fourier_components = 0.01,
           frequency = 2,
           resample = 1,
           time_min = 0,
           time_max = Inf) {
    df_fourier <- df |>
      ## in case we have negative intensity
      ## 100 to be on the safe side
      dplyr::mutate(intensity = intensity + 100) |>
      dplyr::mutate(intensity = intensity - (min(intensity) - 0.001)) |>
      dplyr::mutate(
        intensity_fourier = filter_fft(
          x = intensity,
          components = fourier_components
        )
      )

    f <- approxfun(
      x = df_fourier$time,
      y = df_fourier$intensity_fourier
    )

    time <- seq(
      from = time_min,
      to = min(max(df_fourier$time), time_max),
      by = 1 / (frequency * 60 * resample)
    )

    intensity <- f(seq(
      from = time_min,
      to = min(max(df_fourier$time), time_max),
      by = 1 / (frequency * 60 * resample)
    ))

    intensity_sharpened <- signal_sharpening(time = time, intensity = intensity)

    df_sharpened <-
      data.frame(
        "time" = time[5:length(time)],
        "intensity" = intensity_sharpened
      )

    return(df_sharpened)
  }
