source(file = "R/baseline_chromatogram.R")
source(file = "R/change_intensity_name.R")
source(file = "R/improve_signals_progress.R")
#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
preprocess_chromatograms <- function(detector = "cad",
                                     list = chromatograms_all[c(FALSE, FALSE, TRUE)],
                                     signal_name = "UV.1_CAD_1_0",
                                     shift = CAD_SHIFT) {
  log_debug(x = "preprocessing", detector, "chromatograms")
  log_debug(x = "harmonizing names")
  chromatograms_original <-
    lapply(list, change_intensity_name, signal_name)

  log_debug(x = "improving chromatograms")
  chromatograms_improved <-
    improve_signals_progress(chromatograms_original)

  names(chromatograms_original) <- names
  names(chromatograms_improved) <- names

  chromatograms_original_long <-
    dplyr::bind_rows(chromatograms_original, .id = "id") |>
    dplyr::mutate(time = time + shift) |>
    dplyr::mutate(intensity = intensity - (min(intensity))) |>
    # dplyr::mutate(intensity = intensity / max(intensity)) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()
  chromatograms_improved_long <-
    dplyr::bind_rows(chromatograms_improved, .id = "id") |>
    dplyr::mutate(time = time + shift) |>
    # dplyr::mutate(intensity = intensity / max(intensity)) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()

  # log_debug(x = "plotting improved chromatograms ...")
  # plot_improved <-
  #   plot_chromatogram(df = chromatograms_improved_long, text = detector)

  log_debug(x = "baselining chromatograms")
  chromatograms_baselined <- chromatograms_improved |>
    lapply(FUN = baseline_chromatogram)

  chromatograms_baselined_long <-
    dplyr::bind_rows(chromatograms_baselined, .id = "id") |>
    dplyr::mutate(intensity = intensity - (min(intensity))) |>
    # dplyr::mutate(intensity = intensity / max(intensity)) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()

  # log_debug(x = "plotting baselined chromatograms")
  # plot_baselined <-
  #   plot_chromatogram(df = chromatograms_baselined_long, text = detector)

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
