source(file = "R/baseline_chromatogram.R")
source(file = "R/change_intensity_name.R")
source(file = "R/improve_signals_progress.R")
#' Preprocess chromatograms
#'
#' @include baseline_chromatogram.R
#' @include change_intensity_name.R
#' @include improve_signals_progress.R
#'
#' @param detector Detector
#' @param list List
#' @param signal_name Signal name
#' @param shift Shift
#'
#' @return A list of preprocessed chromatograms
#'
#' @examples NULL
preprocess_chromatograms <- function(detector = "cad",
                                     list = chromatograms_all[c(FALSE, FALSE, TRUE)],
                                     signal_name = "UV.1_CAD_1_0",
                                     shift = CAD_SHIFT) {
  message("preprocessing", detector, "chromatograms")
  message("harmonizing names")
  chromatograms_original <-
    lapply(list, change_intensity_name, signal_name)

  message("improving chromatograms")
  chromatograms_improved <-
    improve_signals_progress(chromatograms_original)

  names(chromatograms_original) <- names
  names(chromatograms_improved) <- names

  chromatograms_original_long <-
    dplyr::bind_rows(chromatograms_original, .id = "id") |>
    dplyr::mutate(time = time + shift) |>
    dplyr::mutate(intensity = intensity - (min(intensity))) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()

  chromatograms_improved_long <-
    dplyr::bind_rows(chromatograms_improved, .id = "id") |>
    dplyr::mutate(time = time + shift) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()

  message("baselining chromatograms")
  chromatograms_baselined <- chromatograms_improved |>
    lapply(FUN = baseline_chromatogram)

  chromatograms_baselined_long <-
    dplyr::bind_rows(chromatograms_baselined, .id = "id") |>
    dplyr::mutate(intensity = intensity - (min(intensity))) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()

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
