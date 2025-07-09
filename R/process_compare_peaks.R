#' Process compare peaks
#'
#' @export
#'
#' @include check_export_dir.R
#' @include compare_peaks.R
#' @include extract_ms_progress.R
#' @include extract_ms_peak.R
#' @include load_chromatograms.R
#' @include load_features.R
#' @include load_name.R
#' @include prepare_features.R
#' @include preprocess_chromatograms.R
#' @include preprocess_peaks.R
#' @include transform_ms.R
#'
#' @param file File path
#' @param features Features path
#' @param type Type. "original", "baselined" or "improved"
#' @param detector Detector
#' @param headers Headers
#' @param export_dir Export directory
#' @param show_example Show example? Default to FALSE
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param min_area Min area
#' @param min_intensity Min intensity
#' @param resample Resample
#' @param shift Shift
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples NULL
process_compare_peaks <- function(
    file = NULL,
    features = NULL,
    type = "baselined",
    detector = "cad",
    headers = c(
      "bpi" = "BasePeak_0",
      "pda" = "PDA#1_TotalAbsorbance_0",
      "cad" = "UV#1_CAD_1_0"
    ),
    export_dir = "data/interim/peaks",
    show_example = FALSE,
    fourier_components = 0.01,
    frequency = 1,
    min_area = 0.005,
    min_intensity = 1E4,
    resample = 1,
    shift = 0.05,
    time_min = 0.5,
    time_max = 32.5) {
  message("loading MS data")
  ms_data <- file |>
    load_ms_data(show_example = show_example)

  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(show_example = show_example, headers = headers)

  message("loading name")
  name <- file |>
    load_name(show_example = show_example)

  message("loading feature table")
  feature_table <- features |>
    load_features(show_example = show_example)

  message("preparing features")
  if (show_example) {
    message("selecting 10 random features for the example")
    set.seed(42)
    feature_table <- feature_table |>
      tidytable::slice_sample(n = 10)
  }
  df_features <- feature_table |>
    prepare_features(min_intensity = min_intensity, name = name)

  message("preprocessing chromatograms")
  switch <- switch(detector,
    "bpi" = headers["bpi"],
    "cad" = headers["cad"],
    "pda" = headers["pda"]
  )
  list <- chromatograms_all[switch |> names()]
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = list,
    # signal_name = signal_name,
    shift = shift,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample
  )

  message("preprocessing peaks")
  peaks_prelist <- preprocess_peaks(
    detector = detector,
    df_features = df_features,
    df_long = switch(type,
      "original" = chromatograms_list$chromatograms_original_long,
      "improved" = chromatograms_list$chromatograms_improved_long,
      "baselined" = chromatograms_list$chromatograms_baselined_long
    ) |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    df_xy = switch(type,
      "original" = chromatograms_list$chromatograms_original[[1]],
      "improved" = chromatograms_list$chromatograms_improved[[1]],
      "baselined" = chromatograms_list$chromatograms_baselined[[1]]
    ),
    min_area = min_area,
    shift = shift,
    name = name
  )

  message("processing ", detector, " peaks")
  message("extracting ms chromatograms (longest step)")
  message(
    "count approx 1 minute per worker per 1000 features (increasing with features number)"
  )
  message("varies a lot depending on features distribution")
  list_ms_chromatograms <- seq_along(
    peaks_prelist$list_df_features_with_peaks_long
  ) |>
    extract_ms_progress(
      ms_data = ms_data,
      rts = peaks_prelist$list_rtr,
      mzs = peaks_prelist$list_mzr,
      nrows = peaks_prelist$list_df_features_with_peaks_long |>
        purrr::map(.f = nrow)
    )

  message("extracting ms peaks")
  list_ms_peaks <- list_ms_chromatograms |>
    purrr::map(
      .f = extract_ms_peak
    )

  message("comparing peaks")
  list_comparison_score <- seq_along(list_ms_peaks) |>
    purrr::map(
      .f = compare_peaks,
      list_ms_peaks = list_ms_peaks,
      peaks_prelist = peaks_prelist
    )

  message("summarizing comparison scores")
  comparison_scores <- list_comparison_score |>
    purrr::flatten()

  message("there are ", length(comparison_scores), " scores calculated")

  message("selecting features with peaks")
  list_df_features_with_scores <- Map(
    function(df, score) {
      df$comparison_score <- score
      df
    },
    peaks_prelist$list_df_features_with_peaks_long,
    comparison_scores
  )

  df_features_with_peaks <- list_df_features_with_scores |>
    tidytable::bind_rows()

  message("final aesthetics")
  df_features_with_peaks_scored <- df_features_with_peaks |>
    tidytable::select(
      sample = id,
      peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = area,
      comparison_score
    ) |>
    tidytable::distinct()

  df_features_without_peaks_scored <-
    peaks_prelist$df_features_without_peaks |>
    tidytable::mutate(comparison_score = NA) |>
    tidytable::select(
      sample = id,
      peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = area,
      comparison_score
    ) |>
    tidytable::distinct()

  message("checking export directory")
  check_export_dir(export_dir)

  message("exporting")

  df_features_with_peaks_scored |>
    tidytable::fwrite(
      file = file.path(
        export_dir,
        name |>
          gsub(pattern = "\\.[^.]+$", replacement = "") |>
          paste("featuresInformed", detector, sep = "_") |>
          paste0(".tsv")
      ),
      sep = "\t"
    )

  df_features_without_peaks_scored |>
    tidytable::fwrite(
      file = file.path(
        export_dir,
        name |>
          gsub(pattern = "\\.[^.]+$", replacement = "") |>
          paste("featuresNotInformed", detector, sep = "_") |>
          paste0(".tsv")
      ),
      sep = "\t"
    )
}
