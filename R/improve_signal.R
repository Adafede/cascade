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
#' @param intensity_offset Offset to add to intensity values to handle negative
#'   intensities. Default is 100. Set to a value larger than the absolute value
#'   of the most negative intensity expected.
#' @param intensity_floor Small value subtracted from minimum intensity to
#'   ensure all values are positive. Default is 0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#'
#' @return A dataframe with improved signal
#'
#' @examples NULL
improve_signal <-
  function(
    df,
    fourier_components = 0.01,
    frequency = 2,
    resample = 1,
    time_min = 0,
    time_max = Inf,
    intensity_offset = 100,
    intensity_floor = 0.001,
    k2 = 250,
    k4 = 1250000,
    sigma = 0.05,
    smoothing_width = 8
  ) {
    df_fourier <- df |>
      ## in case we have negative intensity
      ## intensity_offset to be on the safe side
      tidytable::mutate(intensity = intensity + intensity_offset) |>
      tidytable::mutate(intensity = intensity - (min(intensity) - intensity_floor)) |>
      tidytable::mutate(
        intensity_fourier = filter_fft(
          x = intensity,
          components = fourier_components
        )
      )

    f <- stats::approxfun(
      x = df_fourier$rtime,
      y = df_fourier$intensity_fourier
    )

    time <- seq(
      from = time_min,
      to = min(max(df_fourier$rtime), time_max),
      by = 1 / (frequency * 60 * resample)
    )

    intensity <- f(seq(
      from = time_min,
      to = min(max(df_fourier$rtime), time_max),
      by = 1 / (frequency * 60 * resample)
    ))

    intensity_sharpened <- signal_sharpening(
      time = time,
      intensity = intensity,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      Smoothing_width = smoothing_width
    )

    ## The signal sharpening function removes the first 4 points due to
    ## derivative calculations and smoothing operations
    trim_start <- 5
    return(data.frame(
      "rtime" = time[trim_start:length(time)],
      "intensity" = intensity_sharpened
    ))
  }
