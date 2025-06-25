#' Extract chromatogram
#'
#' @include change_intensity_name.R
#'
#' @param list List
#' @param type Type
#' @param headers Headers
#'
#' @return An extracted chromatogram
#'
#' @examples NULL
extract_chromatogram <- function(list, type, headers) {
  stopifnot(
    "type must be one of 'bpi, 'cad',or 'pda'" = type %in%
      c("bpi", "cad", "pda")
  )
  name <- switch(type,
    "bpi" = headers["bpi"] |> as.character(),
    "cad" = headers["cad"] |> as.character(),
    "pda" = headers["pda"] |> as.character()
  )
  return(
    list[which(headers == name)] |>
      purrr::pluck(1) |>
      change_intensity_name()
  )
}
