#' Signal sharpening
#'
#' @include second_der.R
#'
#' @param time time
#' @param intensity intensity
#' @param k2 K2
#' @param k4 K4
#' @param sigma Sigma
#' @param Smoothing_width Smoothing width
#' @param Baseline_adjust Baseline adjust
#'
#' @return A sharpened signal
#'
#' @examples NULL
signal_sharpening <- function(time,
                              intensity,
                              k2 = 250,
                              k4 = 1250000,
                              sigma = 0.05,
                              Smoothing_width = 8,
                              Baseline_adjust = 0) {
  smooth_1 <- caTools::runmean(
    x = intensity,
    k = Smoothing_width,
    align = "center"
  ) + Baseline_adjust

  smooth_2 <- caTools::runmean(
    x = smooth_1,
    k = Smoothing_width,
    align = "center"
  )

  deriv_2 <- second_der(
    x = time,
    y = smooth_2
  )

  smooth_3 <- caTools::runmean(
    x = deriv_2,
    k = Smoothing_width,
    align = "center"
  )

  deriv_4 <- second_der(
    x = time[3:length(time)],
    y = smooth_3
  )

  smooth_4 <- caTools::runmean(
    x = deriv_4,
    k = 8,
    align = "center"
  )

  sharpened <- smooth_1[5:length(smooth_1)] - (sigma / k2 * smooth_3[3:length(smooth_3)]) + (sigma / k4 * smooth_4)
  sharpened[is.na(sharpened)] <- 0
  # sharpened <- sharpened / max(sharpened)

  return(sharpened)
}
