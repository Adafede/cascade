#' Signal sharpening
#'
#' @include second_der.R
#'
#' @param time time
#' @param intensity intensity
#' @param k2 K2 parameter controlling the weight of the second derivative in
#'   signal sharpening. Default is 250. Lower values increase the sharpening
#'   effect from the second derivative.
#' @param k4 K4 parameter controlling the weight of the fourth derivative in
#'   signal sharpening. Default is 1250000. Lower values increase the
#'   sharpening effect from the fourth derivative.
#' @param sigma Sigma parameter for derivative weighting. Default is 0.05.
#'   Higher values increase the overall sharpening effect.
#' @param Smoothing_width Smoothing width for the running mean filter. Default
#'   is 8. Higher values provide more smoothing but reduce resolution.
#' @param Baseline_adjust Baseline adjustment value. Default is 0.
#'
#' @return A sharpened signal
#'
#' @examples NULL
signal_sharpening <- function(
  time,
  intensity,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  Smoothing_width = 8,
  Baseline_adjust = 0
) {
  smooth_1 <- caTools::runmean(
    x = intensity,
    k = Smoothing_width,
    align = "center"
  ) +
    Baseline_adjust

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
    k = Smoothing_width,
    align = "center"
  )

  sharpened <- smooth_1[5:length(smooth_1)] -
    (sigma / k2 * smooth_3[3:length(smooth_3)]) +
    (sigma / k4 * smooth_4)
  sharpened[is.na(sharpened)] <- 0
  # sharpened <- sharpened / max(sharpened)

  return(sharpened)
}
