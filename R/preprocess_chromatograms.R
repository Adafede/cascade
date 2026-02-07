#' Preprocess chromatograms
#'
#' @include baseline_chromatogram.R
#' @include change_intensity_name.R
#' @include improve_signals_progress.R
#'
#' @param detector Detector
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param list List
#' @param name Name
#' @param resample Resample
#' @param shift Shift
#' @param time_min Time min
#' @param time_max Time max
#' @param intensity_offset Offset to add to intensity values to handle negative
#'   intensities. Default is 100.
#' @param intensity_floor Small value subtracted from minimum intensity.
#'   Default is 0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#' @param baseline_method Method for baseline correction. Default is
#'   "peakDetection". See \code{\link[baseline]{baseline}} for available
#'   methods.
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
  # signal_name = "UV.1_CAD_1_0",
  time_min = 0,
  time_max = Inf,
  intensity_offset = 100,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection"
) {
  message("preprocessing ", detector, " chromatograms")
  # message("harmonizing names")
  # chromatograms_original <-
  #   purrr::map(.x = list, .f = change_intensity_name, name_intensity = signal_name)
  chromatograms_original <- list

  message("improving chromatograms")
  chromatograms_improved <-
    improve_signals_progress(
      xs = chromatograms_original,
      fourier_components = fourier_components,
      frequency = frequency,
      resample = resample,
      time_min = time_min,
      time_max = time_max,
      intensity_offset = intensity_offset,
      intensity_floor = intensity_floor,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      smoothing_width = smoothing_width
    )

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
