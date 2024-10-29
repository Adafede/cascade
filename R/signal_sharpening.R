source(file = "R/second_der.R")

#' Signal sharpening
#'
#' @include second_der.R
#'
#' @param Time Time
#' @param Intensity Intensity
#' @param k2 K2
#' @param k4 K4
#' @param Smoothing_width Smoothing width
#' @param Baseline_adjust Baseline adjust
#'
#' @return A sharpened signal
#'
#' @examples NULL
signal_sharpening <- function(Time = timeow,
                              Intensity = intensityeah,
                              k2 = 250,
                              k4 = 1250000,
                              sigma = 0.05,
                              Smoothing_width = 8,
                              Baseline_adjust = 0) {
  smooth_1 <- zoo::rollmean(
    x = Intensity,
    k = Smoothing_width,
    align = "center",
    fill = 0
  ) + Baseline_adjust

  smooth_2 <- zoo::rollmean(
    x = smooth_1,
    k = Smoothing_width,
    align = "center",
    fill = 0
  )

  deriv_2 <- second_der(
    x = Time,
    y = smooth_2
  )

  smooth_3 <- zoo::rollmean(
    x = deriv_2,
    k = Smoothing_width,
    align = "center",
    fill = 0
  )

  deriv_4 <- second_der(
    x = Time[3:length(Time)],
    y = smooth_3
  )

  smooth_4 <- zoo::rollmean(
    x = deriv_4,
    k = 8,
    align = "center",
    fill = 0
  )

  sharpened <- smooth_1[5:length(smooth_1)] - (sigma / k2 * smooth_3[3:length(smooth_3)]) + (sigma / k4 * smooth_4)
  sharpened[is.na(sharpened)] <- 0
  # sharpened <- sharpened / max(sharpened)

  return(sharpened)
}
