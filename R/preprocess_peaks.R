#' Preprocess peaks
#'
#' @include join_peaks.R
#' @include normalize_chromato.R
#' @include peaks_progress.R
#' @include prepare_mz.R
#' @include prepare_peaks.R
#' @include prepare_rt.R
#'
#' @param detector Detector
#' @param df_features DF features
#' @param df_long DF long
#' @param list List
#' @param name Name
#' @param shift shift
#' @param min_area Minium area
#'
#' @return A list of lists and dataframe with preprocessed peaks
#'
#' @examples NULL
preprocess_peaks <- function(detector = "cad",
                             df_features,
                             df_long,
                             list,
                             name,
                             shift = 0,
                             min_area = 0) {
  message("preprocessing ", detector, " peaks")
  ## data.table call outside of future because buggy else
  peaks <- peaks_progress(list)

  names(peaks) <- name

  ## data.table call outside of future because buggy else
  peaks_long <- dplyr::bind_rows(peaks, .id = "id") |>
    tidytable::data.table()

  message("joining peaks ...")
  df_peaks <-
    join_peaks(
      chromatograms = df_long,
      peaks = peaks_long,
      min_area = min_area
    )

  data.table::setkey(df_peaks, rt_min, rt_max)

  message("joining within given rt tolerance")
  df_features_peaks <-
    data.table::foverlaps(df_features, df_peaks)

  df_features_with_peaks <- df_features_peaks |>
    tidytable::select(-rt_1, -rt_2) |>
    tidytable::filter(!is.na(peak_id)) |>
    tidytable::distinct()

  # df_new_with_cad <- df_new_with_cad |> ## TODO DONT FORGET
  #   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) |> ## TODO DONT FORGET
  #   sample_n(500) ## TODO DONT FORGET

  message("selecting features outside peaks")
  df_features_without_peaks <- df_features_peaks |>
    tidytable::filter(is.na(peak_id)) |>
    tidytable::distinct() # |>
  # filter(sample %in% df_features_with_peaks$sample)

  # df_new_without <- df_new_without |> ## TODO DONT FORGET
  #   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) ## TODO DONT FORGET

  message("splitting by file")
  list_df_features_with_peaks <- df_features_with_peaks |>
    tidytable::group_split(id)

  names(list_df_features_with_peaks) <-
    unique(df_features_with_peaks$id)

  message("splitting by peak")
  list_df_features_with_peaks_per_peak <-
    future.apply::future_lapply(
      X = list_df_features_with_peaks,
      FUN = function(x) {
        x <- x |>
          tidytable::group_split(peak_id)
        return(x)
      }
    )

  list_df_features_with_peaks_long <-
    list_df_features_with_peaks_per_peak |>
    purrr::flatten()

  # message( "retrieving ms files")
  # list_dda_with_peak <-
  #   future_lapply(
  #     X = list_df_features_with_peaks_long,
  #     FUN = function(x) {
  #       x <- x |>
  #         filter_ms(shift = CAD_SHIFT)
  #       return(x)
  #     }
  #   )
  #
  message("normalizing chromato")
  list_chromato_with_peak <-
    future.apply::future_lapply(X = list_df_features_with_peaks_long, FUN = normalize_chromato, list = list)

  message("preparing peaks chromato")
  list_chromato_peaks <-
    future.apply::future_lapply(X = list_chromato_with_peak, FUN = prepare_peaks)

  message("preparing rt")
  list_rtr <-
    future.apply::future_lapply(X = list_df_features_with_peaks_long, FUN = prepare_rt, shift = shift)

  message("preparing mz")
  list_mzr <-
    future.apply::future_lapply(X = list_df_features_with_peaks_long, FUN = prepare_mz)

  returned_list <- list(
    list_df_features_with_peaks_long,
    # list_dda_with_peak,
    list_chromato_peaks,
    list_rtr,
    list_mzr,
    df_features_without_peaks
  )
  names(returned_list) <- c(
    "list_df_features_with_peaks_long",
    # "list_dda_with_peak",
    "list_chromato_peaks",
    "list_rtr",
    "list_mzr",
    "df_features_without_peaks"
  )
  return(returned_list)
}
