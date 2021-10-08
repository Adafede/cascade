source(file = "R/second_der.R")

#' Title
#'
#' @param Time
#' @param Intensity
#' @param K2
#' @param K4
#' @param Smoothing_width
#' @param Baseline_adjust
#'
#' @return
#' @export
#'
#' @examples
signal_sharpening <- function(Time = timeow,
                              Intensity = intensityeah,
                              K2 = k2,
                              K4 = k4,
                              Smoothing_width = smoothing_width,
                              Baseline_adjust = baseline_adjust) {
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
    k = smoothing_width,
    align = "center",
    fill = 0
  )

  sharpened <- smooth_1[5:length(smooth_1)] - (k2 * smooth_3[3:length(smooth_3)]) + (k4 * smooth_4)
  sharpened[is.na(sharpened)] <- 0
  sharpened <- sharpened / max(sharpened)

  return(sharpened)
}
