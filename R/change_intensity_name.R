#' Title
#'
#' @param df
#' @param name
#'
#' @return
#' @export
#'
#' @examples
change_intensity_name <- function(df, name) {
  df |>
    dplyr::select(time,
      intensity = !!as.name(name)
    )
}
