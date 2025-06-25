#' Change intensity name
#'
#' @param df Dataframe
#' @param name_rt Name RT
#' @param name_intensity Name intensity
#'
#' @return A dataframe with changed intensity name
#'
#' @examples NULL
change_intensity_name <- function(
    df,
    name_rt = "rtime",
    name_intensity = "intensity") {
  df |>
    ## Old, see https://github.com/sneumann/mzR/issues/304
    ## TODO Fix could be better
    # tidytable::select(time, intensity = !!as.name(name)) |>
    tidytable::select(
      rtime = !!as.name(name_rt),
      intensity = !!as.name(name_intensity)
    ) |>
    data.frame()
}
