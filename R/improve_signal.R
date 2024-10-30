source(file = "R/signal_sharpening.R")

#' Improve signal
#'
#' @include signal_sharpening.R
#'
#' @param df Dataframe
#' @param fourier_components Fourier components
#' @param time_min Time min
#' @param time_max Time max
#' @param frequency Frequency
#' @param resample Resample
#'
#' @return A dataframe with improved signal
#'
#' @examples NULL
improve_signal <-
  function(df,
           fourier_components = 0.01,
           time_min = 0,
           time_max = Inf,
           frequency = 2,
           resample = 1) {
    df_fourier <- df |>
      ## in case we have negative intensity
      ## 100 to be on the safe side
      dplyr::mutate(intensity = intensity + 100) |>
      dplyr::mutate(intensity = intensity - (min(intensity) - 0.001)) |>
      dplyr::mutate(
        intensity_fourier = nucleR::filterFFT(
          intensity,
          pcKeepComp = fourier_components,
          # pcKeepComp = 0.01,
          showPowerSpec = FALSE,
          useOptim = TRUE
        )
      )

    f <- approxfun(
      x = df_fourier$time,
      y = df_fourier$intensity_fourier
    )

    timeow <<- seq(
      from = time_min,
      to = min(max(df_fourier$time), time_max),
      by = 1 / (frequency * 60 * resample)
    )

    intensityeah <<- f(seq(
      from = time_min,
      to = min(max(df_fourier$time), time_max),
      by = 1 / (frequency * 60 * resample)
    ))

    intensity_sharpened <- signal_sharpening()

    df_sharpened <-
      data.frame(
        "time" = timeow[5:length(timeow)],
        "intensity" = intensity_sharpened
      )

    return(df_sharpened)
  }
