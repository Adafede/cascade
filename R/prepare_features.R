#' Prepare features
#'
#' @param df
#'
#' @return A dataframe of prepared features
#'
#' @export
#'
#' @examples NULL
prepare_features <- function(df, min_intensity) {
  message("... preparing features")
  df <- df |>
    tidytable::select(
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

  message("... keeping features above desired intensity only")
  df_features <- df |>
    tidytable::filter(intensity_max >= min_intensity) |>
    tidytable::data.table()

  message("setting joining keys")
  data.table::setkey(df_features, rt_1, rt_2)

  return(df_features)
}
