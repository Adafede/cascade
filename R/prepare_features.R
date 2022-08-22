#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
prepare_features <- function(df) {
  colnames(df) <-
    gsub(
      pattern = ".Peak.area",
      replacement = "",
      x = colnames(df)
    )

  log_debug(x = "... removing \"row m/z\" and from \"row retention time\" columns")
  df <- df |>
    dplyr::select(
      rt = "row retention time",
      mz = "row m/z",
      everything(),
      # TODO improve if Sirius
      -"correlation group ID",
      -"annotation network number",
      -"best ion",
      -"auto MS2 verify",
      -"identified by n=",
      -"partners",
      -"neutral M mass"
    ) |>
    dplyr::select(-(ncol(df) - 7)) |>
    tibble::column_to_rownames(var = "row ID")

  log_debug(x = "... keeping features above desired intensity only")
  df_features <- df |>
    tibble::rownames_to_column() |>
    dplyr::group_by(rowname, rt, mz) |>
    tidyr::gather(column, value, -rowname, -rt, -mz) |>
    dplyr::ungroup() |>
    dplyr::mutate(column = gsub(
      pattern = "^X",
      replacement = "",
      x = column
    )) |>
    dplyr::arrange(rowname, dplyr::desc(value)) |>
    dplyr::filter(value >= INTENSITY_MS_MIN) |>
    dplyr::mutate(column = gsub(
      pattern = ".Peak.area",
      replacement = "",
      x = column
    )) |>
    dplyr::select(
      feature_id = rowname,
      rt,
      mz,
      sample = column,
      intensity = value
    ) |>
    dplyr::mutate(
      rt_1 = as.numeric(rt),
      rt_2 = as.numeric(rt),
      mz_min = (1 - (1E-6 * PPM)) * as.numeric(mz),
      mz_max = (1 + (1E-6 * PPM)) * as.numeric(mz)
    ) |>
    data.table::data.table()

  log_debug(x = "setting joining keys")
  data.table::setkey(df_features, rt_1, rt_2)

  return(df_features)
}
