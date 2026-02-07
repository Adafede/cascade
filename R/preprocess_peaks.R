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
#' @param df_xy DF X Y
#' @param name Name
#' @param shift shift
#' @param min_area Minimum area
#' @param sd_max Maximum standard deviation for peak filtering. Default is 50.
#' @param max_iter Maximum iterations for peak fitting. Default is 1000.
#' @param noise_threshold Noise threshold for peak detection. Default is 0.001.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw". Default
#'   is "egh".
#' @param intensity_threshold Minimum normalized intensity threshold for
#'   filtering in normalize_chromato. Default is 0.1.
#'
#' @return A list of lists and dataframe with preprocessed peaks
#'
#' @examples NULL
preprocess_peaks <- function(
  detector = "cad",
  df_features,
  df_long,
  df_xy,
  name,
  shift = 0,
  min_area = 0,
  sd_max = 50,
  max_iter = 1000,
  noise_threshold = 0.001,
  fit = "egh",
  intensity_threshold = 0.1
) {
  message("preprocessing ", detector, " peaks")
  ## data.table call outside of future because buggy else
  peaks <- peaks_progress(
    df_xy = df_xy,
    sd_max = sd_max,
    max_iter = max_iter,
    noise_threshold = noise_threshold,
    fit = fit
  )

  ## data.table call outside of future because buggy else
  peaks_long <- tidytable::bind_rows(peaks, .id = "id") |>
    tidytable::data.table()

  message("joining peaks")
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

  message("selecting features outside peaks")
  df_features_without_peaks <- df_features_peaks |>
    tidytable::filter(is.na(peak_id)) |>
    tidytable::distinct()

  message("splitting by file")
  list_df_features_with_peaks <- df_features_with_peaks |>
    tidytable::group_split(id)

  names(list_df_features_with_peaks) <-
    unique(df_features_with_peaks$id)

  message("splitting by peak")
  list_df_features_with_peaks_per_peak <- list_df_features_with_peaks |>
    purrr::map(
      .f = function(x) {
        x <- x |>
          tidytable::group_split(peak_id)
        return(x)
      }
    )

  list_df_features_with_peaks_long <-
    list_df_features_with_peaks_per_peak |>
    purrr::flatten()

  message("normalizing chromato")
  list_chromato_with_peak <- list_df_features_with_peaks_long |>
    purrr::map(
      .f = normalize_chromato,
      df_xy = df_xy,
      intensity_threshold = intensity_threshold
    )

  message("preparing peaks chromato")
  list_chromato_peaks <- list_chromato_with_peak |>
    purrr::map(
      .f = prepare_peaks
    )

  message("preparing rt")
  list_rtr <- list_df_features_with_peaks_long |>
    purrr::map(
      .f = prepare_rt,
      shift = shift
    )

  message("preparing mz")
  list_mzr <- list_df_features_with_peaks_long |>
    purrr::map(
      .f = prepare_mz
    )

  returned_list <- list(
    list_df_features_with_peaks_long,
    list_chromato_peaks,
    list_rtr,
    list_mzr,
    df_features_without_peaks
  )
  names(returned_list) <- c(
    "list_df_features_with_peaks_long",
    "list_chromato_peaks",
    "list_rtr",
    "list_mzr",
    "df_features_without_peaks"
  )
  return(returned_list)
}
