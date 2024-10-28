source(file = "R/join_peaks.R")
source(file = "R/normalize_chromato.R")
source(file = "R/peaks_progress.R")
source(file = "R/prepare_mz.R")
source(file = "R/prepare_peaks.R")
source(file = "R/prepare_rt.R")
#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
preprocess_peaks <- function(detector = "cad",
                             list = chromatograms_list_cad$chromatograms_baselined,
                             df_long = chromatograms_list_cad$chromatograms_baselined_long,
                             area_min) {
  log_debug(x = "preprocessing", detector, "peaks")
  ## data.table call outside of future because buggy else
  peaks <- peaks_progress(list)

  names(peaks) <- names

  ## data.table call outside of future because buggy else
  peaks_long <- dplyr::bind_rows(peaks, .id = "id") |>
    tidytable::data.table()

  log_debug(x = "joining peaks ...")
  df_peaks <-
    join_peaks(chromatograms = df_long, peaks = peaks_long, min_area = area_min)

  data.table::setkey(df_peaks, rt_min, rt_max)

  log_debug(x = "joining within given rt tolerance")
  df_features_peaks <-
    data.table::foverlaps(df_features, df_peaks)

  df_features_with_peaks <- df_features_peaks |>
    tidytable::select(-rt_1, -rt_2) |>
    tidytable::filter(!is.na(peak_id)) |>
    tidytable::distinct()

  # df_new_with_cad <- df_new_with_cad |> ## TODO DONT FORGET
  #   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) |> ## TODO DONT FORGET
  #   sample_n(500) ## TODO DONT FORGET

  log_debug(x = "selecting features outside peaks")
  df_features_without_peaks <- df_features_peaks |>
    tidytable::filter(is.na(peak_id)) |>
    tidytable::distinct() # |>
  # filter(sample %in% df_features_with_peaks$sample)

  # df_new_without <- df_new_without |> ## TODO DONT FORGET
  #   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) ## TODO DONT FORGET

  log_debug(x = "splitting by file")
  list_df_features_with_peaks <- df_features_with_peaks |>
    tidytable::group_split(id)

  names(list_df_features_with_peaks) <-
    unique(df_features_with_peaks$id)

  log_debug(x = "splitting by peak")
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

  # log_debug(x = "retrieving ms files")
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
  log_debug(x = "normalizing chromato")
  list_chromato_with_peak <-
    future.apply::future_lapply(
      X = list_df_features_with_peaks_long,
      FUN = normalize_chromato,
      list = list
    )

  log_debug(x = "preparing peaks chromato")
  list_chromato_peaks <-
    future.apply::future_lapply(
      X = list_chromato_with_peak,
      FUN = prepare_peaks
    )

  log_debug(x = "preparing rt")
  list_rtr <-
    future.apply::future_lapply(
      X = list_df_features_with_peaks_long,
      FUN = prepare_rt
    )

  log_debug(x = "preparing mz")
  list_mzr <-
    future.apply::future_lapply(
      X = list_df_features_with_peaks_long,
      FUN = prepare_mz
    )

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
