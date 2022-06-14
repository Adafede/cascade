#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
preprocess_chromatograms <- function(detector = "cad") {
  list <- switch(detector,
    "bpi" = chromatograms_all[c(TRUE, FALSE, FALSE)],
    "cad" = chromatograms_all[c(FALSE, FALSE, TRUE)],
    "pda" = chromatograms_all[c(FALSE, TRUE, FALSE)]
  )
  signal_name <- switch(detector,
    "bpi" = "BasePeak_0",
    "cad" = "UV.1_CAD_1_0",
    "pda" = "PDA.1_TotalAbsorbance_0"
  )
  shift <- switch(detector,
    "bpi" = 0,
    "cad" = CAD_SHIFT,
    "pda" = PDA_SHIFT
  )

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
    dplyr::mutate(time = time + shift)
  chromatograms_improved_long <-
    dplyr::bind_rows(chromatograms_improved, .id = "id") |>
    dplyr::mutate(time = time + shift)

  # log_debug(x = "plotting improved chromatograms ...")
  # plot_improved <-
  #   plot_chromatogram(df = chromatograms_improved_long, text = detector)

  log_debug(x = "baselining chromatograms ...")
  chromatograms_baselined <-
    transform_baseline(x = chromatograms_improved)

  chromatograms_baselined_long <-
    dplyr::bind_rows(chromatograms_baselined, .id = "id") |>
    dplyr::mutate(intensity = intensity / max(intensity)) |>
    dplyr::mutate(rt_1 = time, rt_2 = time) |>
    data.table::data.table()

  # log_debug(x = "plotting baselined chromatograms ...")
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
