#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
prepare_features <- function(df) {
  log_debug(x = "... preparing features")
  df <- df |>
    dplyr::select(
      feature_id = "id",
      rt = "rt",
      mz = "mz",
      area = "area",
      rt_1 = paste0("datafile:", names, ".mzML:rt_range:min"),
      rt_2 = paste0("datafile:", names, ".mzML:rt_range:max"),
      mz_min = paste0("datafile:", names, ".mzML:mz_range:min"),
      mz_max = paste0("datafile:", names, ".mzML:mz_range:max"),
      intensity_min = paste0("datafile:", names, ".mzML:intensity_range:min"),
      intensity_max = paste0("datafile:", names, ".mzML:intensity_range:max"),
    ) |>
    dplyr::mutate_all(as.numeric)

  log_debug(x = "... keeping features above desired intensity only")
  df_features <- df |>
    dplyr::filter(intensity_max >= INTENSITY_MS_MIN) |>
    data.table::data.table()

  log_debug(x = "setting joining keys")
  data.table::setkey(df_features, rt_1, rt_2)

  return(df_features)
}
