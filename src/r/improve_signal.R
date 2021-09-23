source(file = "src/R/signal_sharpening.R")

#' Title
#'
#' @param df
#' @param fourier_components
#' @param time_min
#' @param time_max
#' @param frequency
#' @param resample
#'
#' @return
#' @export
#'
#' @examples
improve_signal <-
  function(df,
           fourier_components = FOURRIER_COMPONENTS,
           time_min = TIME_MIN,
           time_max = TIME_MAX,
           frequency = FREQUENCY,
           resample = RESAMPLE) {
    df_fourier <- df |>
      dplyr::mutate(intensity = intensity - (min(intensity))) |>
      dplyr::mutate(
        intensity_fourier = nucleR::filterFFT(
          intensity,
          pcKeepComp = fourier_components,
          ## 0.01
          showPowerSpec = FALSE,
          useOptim = TRUE
        )
      )
    
    f <- approxfun(
      x = df_fourier$time,
      y = df_fourier$intensity_fourier
    )
    
    time <- seq(
      from = time_min,
      to = time_max,
      by = 1 / (frequency * 60 * resample)
    )
    
    intensity <- f(seq(
      from = time_min,
      to = time_max,
      by = 1 / (frequency * 60 * resample)
    ))
    
    intensity <- signal_sharpening(Intensity = intensity)
    
    df_sharpened <- data.frame(time, intensity)
    
    return(df_sharpened)
  }