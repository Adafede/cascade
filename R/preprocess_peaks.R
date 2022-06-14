#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
preprocess_peaks <- function(detector = "cad") {
  list <- switch(detector,
    "bpi" = chromatograms_list_bpi$chromatograms_baselined,
    "cad" = chromatograms_list_cad$chromatograms_baselined,
    "pda" = chromatograms_list_pda$chromatograms_baselined
  )
  df_long <- switch(detector,
    "bpi" = chromatograms_list_bpi$chromatograms_baselined_long,
    "cad" = chromatograms_list_cad$chromatograms_baselined_long,
    "pda" = chromatograms_list_pda$chromatograms_baselined_long
  )

  log_debug(x = "preprocessing", detector, "peaks")
  #' data.table call outside of future because buggy else
  peaks <- peaks_progress(list)

  names(peaks) <- names

  #' data.table call outside of future because buggy else
  peaks_long <- dplyr::bind_rows(peaks, .id = "id") |>
    data.table::data.table()

  log_debug(x = "joining peaks ...")
  df_peaks <-
    join_peaks(chromatograms = df_long, peaks = peaks_long)

  data.table::setkey(df_peaks, rt_min, rt_max)

  log_debug(x = "joining within given rt tolerance")
  df_features_peaks <-
    data.table::foverlaps(df_features, df_peaks)

  df_features_with_peaks <- df_features_peaks |>
    dplyr::select(-rt_1, -rt_2) |>
    dplyr::filter(!is.na(peak_id))

  # df_new_with_cad <- df_new_with_cad |> #' TODO DONT FORGET
  #   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) |> #' TODO DONT FORGET
  #   sample_n(500) #' TODO DONT FORGET

  log_debug(x = "selecting features outside peaks")
  df_features_without_peaks <- df_features_peaks |>
    dplyr::filter(is.na(peak_id)) |>
    dplyr::filter(sample %in% df_features_with_peaks$sample)

  # df_new_without <- df_new_without |> #' TODO DONT FORGET
  #   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) #' TODO DONT FORGET

  log_debug(x = "splitting by file")
  list_df_features_with_peaks <- df_features_with_peaks |>
    dplyr::group_split(id)

  names(list_df_features_with_peaks) <-
    unique(df_features_with_peaks$id)

  log_debug(x = "splitting by peak")
  list_df_features_with_peaks_per_peak <-
    parallel::mclapply(
      X = list_df_features_with_peaks,
      FUN = dplyr::group_split,
      peak_id
    )

  list_df_features_with_peaks_long <-
    list_df_features_with_peaks_per_peak |>
    purrr::flatten()

  log_debug(x = "retrieving ms files")
  list_dda_with_peak <-
    parallel::mclapply(
      X = list_df_features_with_peaks_long,
      FUN = filter_ms,
      shift = CAD_SHIFT
    )

  log_debug(x = "normalizing chromato")
  list_chromato_with_peak <-
    parallel::mclapply(
      X = list_df_features_with_peaks_long,
      FUN = normalize_chromato,
      list = list
    )

  log_debug(x = "preparing peaks chromato")
  list_chromato_peaks <-
    parallel::mclapply(
      X = list_chromato_with_peak,
      FUN = prepare_peaks
    )

  log_debug(x = "preparing rt")
  list_rtr <-
    parallel::mclapply(
      X = list_df_features_with_peaks_long,
      FUN = prepare_rt
    )

  log_debug(x = "preparing mz")
  list_mzr <-
    parallel::mclapply(
      X = list_df_features_with_peaks_long,
      FUN = prepare_mz
    )

  returned_list <- list(
    list_df_features_with_peaks_long,
    list_dda_with_peak,
    list_chromato_peaks,
    list_rtr,
    list_mzr,
    df_features_without_peaks
  )
  names(returned_list) <- c(
    "list_df_features_with_peaks_long",
    "list_dda_with_peak",
    "list_chromato_peaks",
    "list_rtr",
    "list_mzr",
    "df_features_without_peaks"
  )
  return(returned_list)
}
