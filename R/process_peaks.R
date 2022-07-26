#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
process_peaks <- function(detector = "cad") {
  log_debug(x = "processing", detector, "peaks")
  log_debug(x = "extracting ms chromatograms (longest step)")
  log_debug(x = "count approx 1 minute per 500 features (increasing with features number)")
  log_debug(x = "varies a lot depending on features distribution")

  list_ms_chromatograms <<-
    parallel::mclapply(
      X = seq_along(switch(detector,
        "bpi" = peaks_prelist_bpi$list_df_features_with_peaks_long,
        "cad" = peaks_prelist_cad$list_df_features_with_peaks_long,
        "pda" = peaks_prelist_pda$list_df_features_with_peaks_long
      )),
      FUN = extract_ms
    )

  log_debug(x = "transforming ms chromatograms")
  list_ms_chromatograms_transformed <<-
    parallel::mclapply(
      X = list_ms_chromatograms,
      FUN = transform_ms
    )

  log_debug(x = "extracting ms peaks")
  list_ms_peaks <<-
    parallel::mclapply(
      X = list_ms_chromatograms_transformed,
      FUN = extract_ms_peak
    )

  log_debug(x = "comparing peaks")
  list_comparison_score <<-
    parallel::mclapply(
      X = seq_along(list_ms_peaks),
      FUN = compare_peaks
    )

  log_debug(x = "selecting features with peaks")
  df_features_with_peaks <<- switch(detector,
    "bpi" = peaks_prelist_bpi$list_df_features_with_peaks_long,
    "cad" = peaks_prelist_cad$list_df_features_with_peaks_long,
    "pda" = peaks_prelist_pda$list_df_features_with_peaks_long
  ) |>
    dplyr::bind_rows()

  log_debug(x = "There are", nrow(df_features_with_peaks), "features with peaks")

  log_debug(x = "summarizing comparison scores")
  comparison_scores <- list_comparison_score |>
    purrr::flatten()

  log_debug(x = "There are", length(comparison_scores), "scores calculated")

  log_debug(x = "joining")
  df_features_with_peaks$comparison_score <-
    as.numeric(comparison_scores)

  log_debug(x = "final aesthetics")
  df_features_with_peaks_scored <<- df_features_with_peaks |>
    dplyr::select(
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
    dplyr::distinct()

  df_features_without_peaks_scored <<- switch(detector,
    "bpi" = peaks_prelist_bpi$df_features_without_peaks,
    "cad" = peaks_prelist_cad$df_features_without_peaks,
    "pda" = peaks_prelist_pda$df_features_without_peaks
  ) |>
    dplyr::mutate(comparison_score = NA) |>
    dplyr::select(
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
    dplyr::distinct()

  returned_list <<- list(
    df_features_with_peaks_scored,
    df_features_without_peaks_scored
  )
  names(returned_list) <- c(
    "df_features_with_peaks_scored",
    "df_features_without_peaks_scored"
  )

  log_debug(x = "checking export directory")
  check_export_dir(EXPORT_DIR)

  log_debug(x = "exporting to ...")
  log_debug(x = EXPORT_DIR)
  readr::write_tsv(
    x = df_features_with_peaks_scored,
    file = file.path(EXPORT_DIR, switch(detector,
      "bpi" = EXPORT_FILE_BPI,
      "cad" = EXPORT_FILE_CAD,
      "pda" = EXPORT_FILE_PDA
    ))
  )
  readr::write_tsv(
    x = df_features_without_peaks_scored,
    file = file.path(EXPORT_DIR, switch(detector,
      "bpi" = EXPORT_FILE_BPI_2,
      "cad" = EXPORT_FILE_CAD_2,
      "pda" = EXPORT_FILE_PDA_2
    ))
  )

  return(returned_list)
}
