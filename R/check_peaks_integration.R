#' Check chromatograms alignment
#'
#' @export
#'
#' @include load_chromatograms.R
#' @include load_features.R
#' @include plot_peak_detection.R
#' @include prepare_features.R
#' @include preprocess_chromatograms.R
#' @include preprocess_peaks.R
#'
#' @param file File path
#' @param features Features path
#' @param chromatogram Chromatogram
#' @param min_area Minimum area
#' @param min_intensity Minimum intensity
#' @param cad_shift CAD shift
#' @param show_example Show example? Default to FALSE
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples
#' \dontrun{
#' check_peaks_integration(show_example = TRUE)
#' }
check_peaks_integration <- function(file = NULL,
                                    features = NULL,
                                    chromatogram = "baselined",
                                    min_area = 0.005,
                                    min_intensity = 1E4,
                                    cad_shift = 0.05,
                                    show_example = FALSE) {
  if (!show_example) {
    name <- file |>
      basename()
  } else {
    name <- c("210619_AR_06_V_03_2_01.mzML")
  }
  message("loading feature table")
  feature_table <- features |>
    load_features(show_example = show_example)

  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(show_example = show_example)

  message("preparing feature list ...")
  df_features <- feature_table |>
    prepare_features(min_intensity = min_intensity, name = name)

  chromatograms_list_cad <- preprocess_chromatograms(
    name = name,
    list = chromatograms_all[c(FALSE, FALSE, TRUE)],
    shift = cad_shift
  )

  peaks <-
    preprocess_peaks(
      df_features = df_features,
      df_long = switch(chromatogram,
        "original" = chromatograms_list_cad$chromatograms_original_long,
        "improved" = chromatograms_list_cad$chromatograms_improved_long,
        "baselined" = chromatograms_list_cad$chromatograms_baselined_long
      ) |>
        tidytable::mutate(intensity = intensity / max(intensity)),
      list = switch(chromatogram,
        "original" = chromatograms_list_cad$chromatograms_original,
        "improved" = chromatograms_list_cad$chromatograms_improved,
        "baselined" = chromatograms_list_cad$chromatograms_baselined
      ),
      min_area = min_area,
      shift = cad_shift,
      name = name
    )

  chromatogram_normalized <- switch(chromatogram,
    "original" = chromatograms_list_cad$chromatograms_original_long |>
      tidytable::bind_rows() |>
      tidytable::filter(row_number() %% 10 == 1) |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    "improved" = chromatograms_list_cad$chromatograms_improved_long |>
      tidytable::bind_rows() |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    "baselined" = chromatograms_list_cad$chromatograms_baselined_long |>
      tidytable::bind_rows() |>
      tidytable::mutate(intensity = intensity / max(intensity))
  )

  peaks_normalized <- peaks$list_df_features_with_peaks_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(
      intensity = intensity_max / max(intensity_max),
      peak_max = peak_max / max(peak_max)
    )

  approx_f <- approxfun(
    x = chromatogram_normalized |>
      tidytable::pull(time),
    y = chromatogram_normalized |>
      tidytable::pull(intensity)
  )

  chromatogram_normalized |>
    plot_peak_detection(df2 = peaks_normalized, fun = approx_f)
}
