#' Process peaks
#'
#' @param detector
#'
#' @return A list of processed peaks
#'
#' @examples NULL
process_peaks <- function(detector = "cad") {
  message("processing", detector, "peaks")
  message("extracting ms chromatograms (longest step)")
  message("count approx 1 minute per 500 features (increasing with features number)")
  message("varies a lot depending on features distribution")

  peaks_prelist <<- switch(detector,
    "bpi" = peaks_prelist_bpi,
    "cad" = peaks_prelist_cad,
    "pda" = peaks_prelist_pda
  )

  list_ms_chromatograms <<-
    extract_ms_progress(xs = seq_along(peaks_prelist$list_df_features_with_peaks_long))

  message("transforming ms chromatograms")
  list_ms_chromatograms_transformed <<-
    future_lapply(
      X = list_ms_chromatograms,
      FUN = transform_ms
    )

  message("extracting ms peaks")
  list_ms_peaks <<-
    future_lapply(
      X = list_ms_chromatograms_transformed,
      FUN = extract_ms_peak
    )

  message("comparing peaks")
  list_comparison_score <<-
    future_lapply(
      X = seq_along(list_ms_peaks),
      FUN = compare_peaks
    )

  message("selecting features with peaks")
  df_features_with_peaks <<-
    peaks_prelist$list_df_features_with_peaks_long |>
    bind_rows()

  message("There are", nrow(df_features_with_peaks), "features with peaks")

  message("summarizing comparison scores")
  comparison_scores <- list_comparison_score |>
    flatten()

  message("There are", length(comparison_scores), "scores calculated")

  message("joining")
  df_features_with_peaks$comparison_score <-
    as.numeric(comparison_scores)

  message("final aesthetics")
  df_features_with_peaks_scored <<- df_features_with_peaks |>
    select(
      sample = id,
      peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = intensity,
      comparison_score
    ) |>
    distinct()

  df_features_without_peaks_scored <<-
    peaks_prelist$df_features_without_peaks |>
    mutate(comparison_score = NA) |>
    select(
      sample = id,
      peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = intensity,
      comparison_score
    ) |>
    distinct()

  returned_list <<- list(
    df_features_with_peaks_scored,
    df_features_without_peaks_scored
  )
  names(returned_list) <- c(
    "df_features_with_peaks_scored",
    "df_features_without_peaks_scored"
  )

  message("checking export directory")
  check_export_dir(EXPORT_DIR)

  message("exporting to ...")
  message(EXPORT_DIR)
  write_tsv(
    x = df_features_with_peaks_scored,
    file = file.path(EXPORT_DIR, switch(detector,
      "bpi" = EXPORT_FILE_BPI,
      "cad" = EXPORT_FILE_CAD,
      "pda" = EXPORT_FILE_PDA
    ))
  )
  write_tsv(
    x = df_features_without_peaks_scored,
    file = file.path(EXPORT_DIR, switch(detector,
      "bpi" = EXPORT_FILE_BPI_2,
      "cad" = EXPORT_FILE_CAD_2,
      "pda" = EXPORT_FILE_PDA_2
    ))
  )

  return(returned_list)
}
