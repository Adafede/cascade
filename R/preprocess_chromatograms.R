#' Preprocess chromatograms
#'
#' @include baseline_chromatogram.R
#' @include change_intensity_name.R
#' @include improve_signals_progress.R
#'
#' @param detector Detector type (e.g., "cad", "bpi", "pda")
#' @param fourier_components Fraction of Fourier components to keep. Default is
#'   0.01.
#' @param frequency Acquisition frequency in Hz. Default is 2.
#' @param list List of chromatograms
#' @param name Sample name(s)
#' @param resample Resampling factor. Default is 1.
#' @param shift Time shift in minutes. Default is 0.
#' @param time_min Time min in minutes. Default is 0.
#' @param time_max Time max in minutes. Default is Inf.
#' @param intensity_floor Small positive value for intensity floor. Default is
#'   0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#' @param baseline_method Method for baseline correction. Default is
#'   "peakDetection". See \code{\link[baseline]{baseline}} for available
#'   methods.
#' @param improve_signal Logical. Whether to apply signal improvement (Fourier
#'   filtering and sharpening). Default is TRUE. Set to FALSE to skip signal
#'   improvement and use original chromatograms.
#'
#' @return A list of preprocessed chromatograms
#'
#' @examples NULL
preprocess_chromatograms <- function(
  detector = "cad",
  fourier_components = 0.01,
  frequency = 2,
  list,
  name,
  resample = 1,
  shift = 0,
  time_min = 0,
  time_max = Inf,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection",
  improve_signal = TRUE
) {
  message("preprocessing ", detector, " chromatograms")
  chromatograms_original <- list

  if (improve_signal) {
    message("improving chromatograms")
    chromatograms_improved <-
      improve_signals_progress(
        xs = chromatograms_original,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        time_min = time_min,
        time_max = time_max,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width
      )
  } else {
    message("skipping signal improvement (using original chromatograms)")
    chromatograms_improved <- chromatograms_original |>
      purrr::map(.f = function(x) {
        x |>
          tidytable::select(rtime, intensity) |>
          tidytable::filter(rtime >= time_min & rtime <= time_max) |>
          data.frame()
      })
  }

  names(chromatograms_original) <- name
  names(chromatograms_improved) <- name

  chromatograms_original_long <-
    tidytable::bind_rows(chromatograms_original, .id = "id") |>
    tidytable::mutate(rtime = rtime + shift) |>
    tidytable::mutate(intensity = intensity - (min(intensity))) |>
    tidytable::mutate(rt_1 = rtime, rt_2 = rtime) |>
    tidytable::data.table()

  chromatograms_improved_long <-
    tidytable::bind_rows(chromatograms_improved, .id = "id") |>
    tidytable::mutate(rtime = rtime + shift) |>
    tidytable::mutate(rt_1 = rtime, rt_2 = rtime) |>
    tidytable::data.table()

  message("baselining chromatograms")
  chromatograms_baselined <- chromatograms_improved |>
    purrr::map(.f = baseline_chromatogram, method = baseline_method)

  chromatograms_baselined_long <-
    tidytable::bind_rows(chromatograms_baselined, .id = "id") |>
    tidytable::mutate(intensity = intensity - (min(intensity))) |>
    tidytable::mutate(rt_1 = rtime, rt_2 = rtime) |>
    tidytable::data.table()

  returned_list <- list(
    chromatograms_original,
    chromatograms_original_long,
    chromatograms_improved,
    chromatograms_improved_long,
    chromatograms_baselined,
    chromatograms_baselined_long
  )
  names(returned_list) <- c(
    "chromatograms_original",
    "chromatograms_original_long",
    "chromatograms_improved",
    "chromatograms_improved_long",
    "chromatograms_baselined",
    "chromatograms_baselined_long"
  )
  return(returned_list)
}
