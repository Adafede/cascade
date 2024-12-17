#' Prepare features
#'
#' @param df Df
#' @param min_intensity Min intensity
#' @param name Name
#'
#' @return A dataframe of prepared features
#'
#' @examples NULL
prepare_features <- function(df, min_intensity, name) {
  message("... preparing features")
  df <- df |>
    tidytable::select(tidytable::all_of(
      c(
        feature_id = "id",
        rt = "rt",
        mz = "mz",
        area = "area",
        status = paste0("datafile:", name, ":feature_state"),
        rt_1 = paste0("datafile:", name, ":rt_range:min"),
        rt_2 = paste0("datafile:", name, ":rt_range:max"),
        mz_min = paste0("datafile:", name, ":mz_range:min"),
        mz_max = paste0("datafile:", name, ":mz_range:max"),
        intensity_min = paste0("datafile:", name, ":intensity_range:min"),
        intensity_max = paste0("datafile:", name, ":intensity_range:max")
      )
    )) |>
    tidytable::filter(status == "DETECTED") |>
    tidytable::select(-status) |>
    tidytable::mutate(tidytable::across(tidytable::everything(), as.numeric))

  message("... keeping features above desired intensity only")
  df_features <- df |>
    tidytable::filter(intensity_max >= min_intensity) |>
    tidytable::data.table()

  message("setting joining keys")
  data.table::setkey(df_features, rt_1, rt_2)

  return(df_features)
}
