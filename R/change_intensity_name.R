#' Change intensity name
#'
#' @param df Dataframe
#' @param name Name
#'
#' @return A dataframe with changed intensity name
#'
#' @examples NULL
change_intensity_name <- function(df, name) {
  df |>
    tidytable::select(time, intensity = !!as.name(name)) |>
    data.frame()
}
