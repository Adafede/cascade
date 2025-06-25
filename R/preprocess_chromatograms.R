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
#' @param signal_name Signal name
#' @param time_min Time min
#' @param time_max Time max
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
  signal_name = "UV.1_CAD_1_0",
  time_min = 0,
  time_max = Inf
) {
  message("preprocessing ", detector, " chromatograms")
  message("harmonizing names")
  chromatograms_original <-
    purrr::map(.x = list, .f = change_intensity_name, name = signal_name)

  message("improving chromatograms")
  chromatograms_improved <-
    improve_signals_progress(
      xs = chromatograms_original,
      fourier_components = fourier_components,
      frequency = frequency,
      resample = resample,
      time_min = time_min,
      time_max = time_max
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
    purrr::map(.f = baseline_chromatogram)

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
